#define _XOPEN_SOURCE
#include <boost/bind.hpp>
#include <graphlab/parallel/fiber.hpp>
//#include <valgrind/valgrind.h>
namespace graphlab {

bool fiber_group::tls_created = false;
pthread_key_t fiber_group::tlskey;

fiber_group::fiber_group(size_t nworkers, size_t stacksize)
    :nworkers(nworkers),
    stacksize(stacksize),
    stop_workers(false),
    flsdeleter(NULL) {
  // initialize the thread local storage keys
  if (!tls_created) {
    pthread_key_create(&tlskey, fiber_group::tls_deleter);
    tls_created = true;
  }

  // set up the queues.
  active_head.next = NULL;
  active_tail = &active_head;
  nactive = 0;
  workers_waiting = 0;
  // launch the workers
  for (size_t i = 0;i < nworkers; ++i) {
    workers.launch(boost::bind(&fiber_group::worker_init, this, i));
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
  if (value->scheduleable) {
    value->next = NULL;
    active_tail->next = value;
    active_tail = value;
    ++nactive;
  }
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

void fiber_group::worker_init(size_t workerid) {
  // create a root context
  create_tls_ptr();
  // set up the tls structure
  tls* t = get_tls_ptr();
  t->prev_fiber = NULL;
  t->cur_fiber = NULL;
  t->garbage = NULL;
  t->worker_id = workerid;
  t->parent = this;
  active_lock.lock();
  while(!stop_workers) {
    fiber* next_fib = t->parent->active_queue_remove();
    if (next_fib != NULL) {
      active_lock.unlock();
      yield_to(next_fib);
      active_lock.lock();
    } else {
      ++workers_waiting;
      active_cond.wait(active_lock);
      --workers_waiting;
    }
  }
  active_lock.unlock();
}

struct trampoline_args {
  boost::function<void(void)> fn;
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
    args->fn();
  } catch (...) {
  }
  delete args;
  fiber_group::exit();
}

size_t fiber_group::launch(boost::function<void(void)> fn) {
  // allocate a stack
  fiber* fib = new fiber;
  fib->parent = this;
  fib->stack = malloc(stacksize);
  fib->id = fiber_id_counter.inc();
  //VALGRIND_STACK_REGISTER(fib->stack, (char*)fib->stack + stacksize);
  fib->fls = NULL;
  fib->next = NULL;
  fib->deschedule_lock = NULL;
  fib->terminate = false;
  fib->descheduled = false;
  fib->scheduleable = true;
  // construct the initial context
  trampoline_args* args = new trampoline_args;
  args->fn = fn;
  fib->initial_trampoline_args = (intptr_t)(args);
  // stack grows downwards.
  fib->context = boost::context::make_fcontext((char*)fib->stack + stacksize,
                                               stacksize,
                                               trampoline);
  fibers_active.inc();

  active_lock.lock();
  active_queue_insert(fib);
  if (workers_waiting) active_cond.signal();
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
    // am I meant to be terminated? or descheduled?
    if (t->cur_fiber &&
        (t->cur_fiber->terminate || t->cur_fiber->descheduled) ) {
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
  fib->lock.lock();
  if (!fib->terminate && !fib->descheduled) {
    fib->lock.unlock();
    // we reschedule it
    // Re-lock the queue
    //printf("Reinserting %ld\n", fib->id);
    active_lock.lock();
    active_queue_insert(fib);
    if (workers_waiting) active_cond.signal();
    active_lock.unlock();
  } else if (fib->descheduled) {
    // unflag descheduled and unset scheduleable
    fib->descheduled = false;
    fib->scheduleable = false;
    if (fib->deschedule_lock) pthread_mutex_unlock(fib->deschedule_lock);
    fib->deschedule_lock = NULL;
    //printf("Descheduling complete %ld\n", fib->id);
    fib->lock.unlock();
  } else if (fib->terminate) {
    fib->lock.unlock();
    // previous fiber is dead. destroy it
    free(fib->stack);
    //VALGRIND_STACK_DEREGISTER(fib->stack);
    // delete the fiber local storage if any
    if (fib->fls && flsdeleter) flsdeleter(fib->fls);
    delete fib;
    // if we are out of threads, signal the join
    if (fibers_active.dec() == 0) {
      join_lock.lock();
      join_cond.signal();
      join_lock.unlock();
    }
  } else {
    // impossible condition
    assert(false);
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

size_t fiber_group::get_tid() {
  return reinterpret_cast<size_t>(get_tls_ptr()->cur_fiber);
}

void fiber_group::deschedule_self(pthread_mutex_t* lock) {
  fiber* fib = get_tls_ptr()->cur_fiber;
  fib->lock.lock();
  assert(fib->descheduled == false);
  assert(fib->scheduleable == true);
  fib->deschedule_lock = lock;
  fib->descheduled = true;
  //printf("Descheduling requested %ld\n", fib->id);
  fib->lock.unlock();
  yield();
}

size_t fiber_group::get_worker_id() {
  fiber_group::tls* tls = get_tls_ptr();
  return tls->worker_id;
}

void fiber_group::schedule_tid(size_t tid) {
  fiber* fib = reinterpret_cast<fiber*>(tid);
  fib->lock.lock();
  // we MUST get here only after the thread was completely descheduled
  // or no deschedule operation has happened yet.
  assert(fib->descheduled == false);
  fib->descheduled = false;
  if (fib->scheduleable == false) {
    // if this thread was descheduled completely. Reschedule it.
    //printf("Scheduling requested %ld\n", fib->id);
    fib->scheduleable = true;
    fib->lock.unlock();
    fib->parent->reschedule_fiber(fib);
  } else {
    //printf("Scheduling requested of running thread %ld\n", fib->id);
    fib->lock.unlock();
  }
}


void fiber_group::set_tls_deleter(void (*deleter)(void*)) {
  flsdeleter = deleter;
}

void* fiber_group::get_tls() {
  fiber_group::tls* f = get_tls_ptr();
  if (f != NULL) {
    return f->cur_fiber->fls;
  } else {
    // cannot get TLS of a non-fiber
    assert(false);
    return NULL;
  }
}

void fiber_group::set_tls(void* tls) {
  fiber_group::tls* f = get_tls_ptr();
  if (f != NULL) {
    f->cur_fiber->fls = tls;
  } else {
    // cannot get TLS of a non-fiber
    assert(false);
  }
}


}
