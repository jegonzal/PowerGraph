#define _XOPEN_SOURCE
#include <boost/bind.hpp>
#include <graphlab/parallel/fiber.hpp>
#include <valgrind/valgrind.h>
namespace graphlab {

bool fiber_group::tls_created = false;
pthread_key_t fiber_group::tlskey;

fiber_group::fiber_group(size_t nworkers, size_t stacksize)
    :nworkers(nworkers),
    stacksize(stacksize),
    stop_workers(false) {
  // initialize the thread local storage keys
  if (!tls_created) {
    pthread_key_create(&tlskey, fiber_group::tls_deleter);
    tls_created = true;
  }

  // set up the queues.
  active_head.next = NULL;
  active_tail = &active_head;
  nactive = 0;

  // launch the workers
  for (size_t i = 0;i < nworkers; ++i) {
    workers.launch(boost::bind(&fiber_group::worker_init, this));
  }
}

fiber_group::~fiber_group() {
  join();
  stop_workers = true;
  active_cond.broadcast();
  workers.join();
  pthread_key_delete(tlskey);
}


void fiber_group::tls_deleter(void* f) {
  fiber_group::tls* t = (fiber_group::tls*)(f);
  delete t;
}

void fiber_group::create_tls_ptr() {
  pthread_setspecific(tlskey, (void*)(new fiber_group::tls));
}


fiber_group::tls* fiber_group::get_tls_ptr() {
  return (fiber_group::tls*) pthread_getspecific(tlskey);
}

fiber_group::fiber* fiber_group::get_active_fiber() {
  return get_tls_ptr()->cur_fiber;
}




void fiber_group::active_queue_insert(fiber_group::fiber* value) {
  value->next = NULL;
  active_tail->next = value;
  active_tail = value;
  ++nactive;
  // might want to handle the signalling mechanism here too
}

fiber_group::fiber* fiber_group::active_queue_remove() {
  fiber_group::fiber* ret = active_head.next;
  if (ret != NULL) {
    active_head.next = ret->next;
    --nactive;
    ret->next = NULL;
    if (active_tail == ret) active_tail = &active_head;
  }
  return ret;
}

void fiber_group::exit() {
  fiber* fib = get_active_fiber();
  if (fib != NULL) {
    // add to garbage.
    fib->terminate = true;
    yield(); // never returns
    assert(false);
  }
}

void fiber_group::worker_init() {
  // create a root context
  create_tls_ptr();
  // set up the tls structure
  tls* t = get_tls_ptr();
  t->prev_fiber = NULL;
  t->cur_fiber = NULL;
  t->garbage = NULL;
  t->parent = this;
  active_lock.lock();
  while(!stop_workers) {
    fiber* next_fib = t->parent->active_queue_remove();
    if (next_fib != NULL) {
      active_lock.unlock();
      yield_to(next_fib);
      active_lock.lock();
    } else {
      active_cond.wait(active_lock);
    }
  }
  active_lock.unlock();
}

struct trampoline_args {
  void (*fn)(void* arg);
  void* param;
};

// the trampoline to call the user function. This function never returns
void fiber_group::trampoline(intptr_t _args) {
  // we may have launched to here by switching in from another fiber.
  // we will need to clean up the previous fiber
  tls* t = get_tls_ptr();
  if (t->prev_fiber) t->parent->reschedule_fiber(t->prev_fiber);
  t->prev_fiber = NULL;

  trampoline_args* args = reinterpret_cast<trampoline_args*>(_args);
  try {
    (args->fn)(args->param);
  } catch (...) {
  }
  delete args;
  fiber_group::exit();
}

size_t fiber_group::launch(void fn(void*), void* param) {
  // allocate a stack
  fiber* fib = new fiber;
  fib->stack = malloc(stacksize);
  fib->id = fiber_id_counter.inc();
  VALGRIND_STACK_REGISTER(fib->stack, (char*)fib->stack + stacksize);
  fib->tls = NULL;
  fib->next = NULL;
  fib->terminate = false;
  // construct the initial context
  trampoline_args* args = new trampoline_args;
  args->fn = fn;
  args->param = param;
  fib->initial_trampoline_args = (intptr_t)(args);
  // stack grows downwards.
  fib->context = boost::context::make_fcontext((char*)fib->stack + stacksize,
                                               stacksize,
                                               trampoline);
  fibers_active.inc();

  active_lock.lock();
  active_queue_insert(fib);
  active_cond.signal();
  active_lock.unlock();
  return reinterpret_cast<size_t>(fib);
}

void fiber_group::yield_to(fiber* next_fib) {
  // the core scheduling logic
  tls* t = get_tls_ptr();
/*
  if (next_fib) {
    printf("yield to: %ld\n", next_fib->id);
    if (t->cur_fiber) {
      printf("from: %ld\n", t->cur_fiber->id);
    }
  } */
  if (next_fib != NULL) {
    // current fiber moves to previous
    // next fiber move to current
    t->prev_fiber = t->cur_fiber;
    t->cur_fiber = next_fib;
    if (t->prev_fiber != NULL) {
      // context switch to fib outside the lock
      boost::context::jump_fcontext(t->prev_fiber->context,
                                    t->cur_fiber->context,
                                    t->cur_fiber->initial_trampoline_args);
    } else {
      boost::context::jump_fcontext(&t->base_context,
                                    t->cur_fiber->context,
                                    t->cur_fiber->initial_trampoline_args);
    }
  } else {
    // ok. there isn't anything to schedule to
    // am I meant to be terminated?
    if (t->cur_fiber && t->cur_fiber->terminate) {
      // yup. killing current fiber
      // context switch back to basecontext which will
      // do the cleanup
      //
      // current fiber moves to previous
      // next fiber (base context) move to current
      // (as identifibed by cur_fiber = NULL)
      t->prev_fiber = t->cur_fiber;
      t->cur_fiber = NULL;
      boost::context::jump_fcontext(t->prev_fiber->context,
                                    &t->base_context,
                                    0);
    } else {
      // nothing to do, and not terminating...
      // then don't yield!
      return;
    }
  }
  // reread the tls pointer because we may have woken up in a different thread
  t = get_tls_ptr();
  if (t->prev_fiber) reschedule_fiber(t->prev_fiber);
  t->prev_fiber = NULL;
}

void fiber_group::reschedule_fiber(fiber* fib) {
  if (!fib->terminate) {
    // we reschedule it
    // Re-lock the queue
    active_lock.lock();
    //printf("Reinserting %ld\n", fib->id);
    active_queue_insert(fib);
    active_cond.signal();
    active_lock.unlock();
  } else {
    // previous fiber is dead. destroy it
    free(fib->stack);
    VALGRIND_STACK_DEREGISTER(fib->stack);
    delete fib;
    // if we are out of threads, signal the join
    if (fibers_active.dec() == 0) {
      join_lock.lock();
      join_cond.signal();
      join_lock.unlock();
    }
  }
}

void fiber_group::yield() {
  // the core scheduling logic
  tls* t = get_tls_ptr();

  // remove some other work to do.
  t->parent->active_lock.lock();
  fiber* next_fib = t->parent->active_queue_remove();
  t->parent->active_lock.unlock();
  t->parent->yield_to(next_fib);
}

void fiber_group::join() {
  join_lock.lock();
  while(fibers_active.value > 0) {
    join_cond.wait(join_lock);
  }
  join_lock.unlock();
}

}
