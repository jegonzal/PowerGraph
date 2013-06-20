#include <boost/bind.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/parallel/fiber_control.hpp>
#include <graphlab/logger/assertions.hpp>
//#include <valgrind/valgrind.h>
namespace graphlab {

bool fiber_control::tls_created = false;
bool fiber_control::instance_created = false;
size_t fiber_control::instance_construct_params_nworkers = 0;
size_t fiber_control::instance_construct_params_affinity_base = 0;
pthread_key_t fiber_control::tlskey;

fiber_control::fiber_control(size_t nworkers, 
                             size_t affinity_base)
    :nworkers(nworkers),
    stop_workers(false),
    flsdeleter(NULL) {
  // initialize the thread local storage keys
  if (!tls_created) {
    pthread_key_create(&tlskey, fiber_control::tls_deleter);
    tls_created = true;
  }

  // set up the queues.
  schedule.resize(nworkers);
  for (size_t i = 0;i < nworkers; ++i) {
    schedule[i].active_head.next = NULL;
    schedule[i].active_tail = &schedule[i].active_head;
    schedule[i].nactive = 0;
  }
  nactive = 0;
  // launch the workers
  for (size_t i = 0;i < nworkers; ++i) {
    workers.launch(boost::bind(&fiber_control::worker_init, this, i), 
                   affinity_base + i);
  }
}

fiber_control::~fiber_control() {
  join();
  stop_workers = true;
  for (size_t i = 0;i < nworkers; ++i) {
    schedule[i].active_lock.lock();
    schedule[i].active_cond.broadcast();
    schedule[i].active_lock.unlock();
  }
  workers.join();

  for (size_t i = 0;i < nworkers; ++i) {
    ASSERT_EQ(schedule[i].nactive, 0);
  }
  ASSERT_EQ(nactive.value, 0);
  pthread_key_delete(tlskey);
}


void fiber_control::tls_deleter(void* f) {
  fiber_control::tls* t = (fiber_control::tls*)(f);
  delete t;
}

void fiber_control::create_tls_ptr() {
  pthread_setspecific(tlskey, (void*)(new fiber_control::tls));
}


fiber_control::tls* fiber_control::get_tls_ptr() {
  return (fiber_control::tls*) pthread_getspecific(tlskey);
}

fiber_control::fiber* fiber_control::get_active_fiber() {
  tls* t = get_tls_ptr();
  if (t != NULL) return t->cur_fiber;
  else return NULL;
}




void fiber_control::active_queue_insert_tail(size_t workerid, fiber_control::fiber* value) {
  if (value->scheduleable) {
    value->next = NULL;
    schedule[workerid].active_tail->next = value;
    schedule[workerid].active_tail = value;
    ++schedule[workerid].nactive;
    ++nactive;
  }
  // might want to handle the signalling mechanism here too
}

void fiber_control::active_queue_insert_head(size_t workerid, fiber_control::fiber* value) {
  if (value->scheduleable) {
    value->next = schedule[workerid].active_head.next;
    schedule[workerid].active_head.next = value;
    // fixup the tail if it was pointing to the head
    if (schedule[workerid].active_tail == &schedule[workerid].active_head) {
      schedule[workerid].active_tail = value;
    }
    ++schedule[workerid].nactive;
    ++nactive;
  }
  // might want to handle the signalling mechanism here too
}

fiber_control::fiber* fiber_control::active_queue_remove(size_t workerid) {
  fiber_control::fiber* ret = schedule[workerid].active_head.next;
  if (ret != NULL) {
    schedule[workerid].active_head.next = ret->next;
    --nactive;
    --schedule[workerid].nactive;
    ret->next = NULL;
    if (schedule[workerid].active_tail == ret) {
      schedule[workerid].active_tail = &schedule[workerid].active_head;
    }
  }
  return ret;
}

void fiber_control::exit() {
  fiber* fib = get_active_fiber();
  if (fib != NULL) {
    // add to garbage.
    fib->terminate = true;
    yield(); // never returns
    ASSERT_MSG(false, "Impossible Condition. Dead Fiber woke up");
  } else {
    ASSERT_MSG(false, "Calling fiber exit not from a fiber");
  }
}

void fiber_control::worker_init(size_t workerid) {
  /*
   * This is the "root" stack for each worker.
   * When there are active user threads associated with this worker, 
   * it will switch directly between the fibers.
   * But, when the worker has no other fiber to run, it will return to this
   * stack and and wait in a condition variable
   */
  // create a root context
  create_tls_ptr();
  // set up the tls structure
  tls* t = get_tls_ptr();
  t->prev_fiber = NULL;
  t->cur_fiber = NULL;
  t->garbage = NULL;
  t->workerid = workerid;
  t->parent = this;

  schedule[workerid].active_lock.lock();
  while(!stop_workers) {
    // get a fiber to run
    fiber* next_fib = t->parent->active_queue_remove(workerid);
    if (next_fib != NULL) {
      // if there is a fiber. yield to it
      schedule[workerid].active_lock.unlock();
      yield_to(next_fib);
      schedule[workerid].active_lock.lock();
    } else {
      // if there is no fiber. wait.
      schedule[workerid].waiting = true;
      schedule[workerid].active_cond.wait(schedule[workerid].active_lock);
      schedule[workerid].waiting = false;
    }
  }
  schedule[workerid].active_lock.unlock();
}

struct trampoline_args {
  boost::function<void(void)> fn;
};

// the trampoline to call the user function. This function never returns
void fiber_control::trampoline(intptr_t _args) {
  // we may have launched to here by switching in from another fiber.
  // we will need to clean up the previous fiber
  tls* t = get_tls_ptr();
  if (t->prev_fiber) t->parent->reschedule_fiber(t->workerid, t->prev_fiber);
  t->prev_fiber = NULL;

  trampoline_args* args = reinterpret_cast<trampoline_args*>(_args);
  try {
    args->fn();
  } catch (...) {
  }
  delete args;
  fiber_control::exit();
}

size_t fiber_control::launch(boost::function<void(void)> fn, 
                             size_t stacksize, 
                             int affinity) {
  // allocate a stack
  fiber* fib = new fiber;
  fib->parent = this;
  fib->stack = malloc(stacksize);
  fib->id = fiber_id_counter.inc();
  fib->affinity = affinity;
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

  // find a place to put the thread
  // pick 2 random numbers. use the choice of 2
  // rb uses a linear congruential generator
  size_t choice = 0;
  if(affinity >= 0) choice = affinity ;
  else choice = get_worker_id();
  if (choice == (size_t)(-1)) {
    choice = load_balanced_worker_choice(fib->id);
  }
  schedule[choice].active_lock.lock();
  active_queue_insert_tail(choice, fib);
  if (schedule[choice].waiting) schedule[choice].active_cond.signal();
  schedule[choice].active_lock.unlock();
  return reinterpret_cast<size_t>(fib);
}

size_t fiber_control::load_balanced_worker_choice(size_t seed) {
  size_t ra = seed % nworkers;
  size_t rb = graphlab::random::fast_uniform<size_t>(0,nworkers-1);
  size_t choice = (schedule[ra].nactive <= schedule[rb].nactive) ? ra : rb;
  return choice;
}

void fiber_control::yield_to(fiber* next_fib) {
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
    // reset the priority flag
    next_fib->priority = false;
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
  // clear the priority flag if set
  if (t->prev_fiber) reschedule_fiber(t->prev_fiber->affinity >= 0 ?
                                         t->prev_fiber->affinity :
                                         t->workerid,
                                      t->prev_fiber);
  t->prev_fiber = NULL;
}

void fiber_control::reschedule_fiber(size_t workerid, fiber* fib) {
  fib->lock.lock();
  if (!fib->terminate && !fib->descheduled) {
    fib->lock.unlock();
    // we reschedule it
    // Re-lock the queue
    //printf("Reinserting %ld\n", fib->id);
    schedule[workerid].active_lock.lock();
    if (!fib->priority) active_queue_insert_tail(workerid, fib);
    else active_queue_insert_head(workerid, fib);
    if (schedule[workerid].waiting) schedule[workerid].active_cond.signal();
    schedule[workerid].active_lock.unlock();
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

void fiber_control::yield() {
  // the core scheduling logic
  tls* t = get_tls_ptr();
  if (t == NULL) return;
  // remove some other work to do.
  fiber_control* parentgroup = t->parent;
  size_t workerid = t->workerid;
  fiber* next_fib = NULL;
  if (parentgroup->schedule[workerid].nactive > 0) {
    parentgroup->schedule[workerid].active_lock.lock();
    next_fib = parentgroup->active_queue_remove(workerid);
    parentgroup->schedule[workerid].active_lock.unlock();
  }
  // no work on my queue!
  if (next_fib == NULL) {
    // ok. do a full sweep. Try to steal some work
    for (size_t i = 1;i < parentgroup->nworkers; ++i) {
      size_t probe = (i + workerid) % parentgroup->nworkers;
      if (parentgroup->schedule[probe].nactive > 0) {
        parentgroup->schedule[probe].active_lock.lock();
        next_fib = parentgroup->active_queue_remove(probe);
        parentgroup->schedule[probe].active_lock.unlock();
        if (next_fib) {
          break;
        }
      }
    }
  }
  t->parent->yield_to(next_fib);
}

void fiber_control::join() {
  join_lock.lock();
  while(fibers_active.value > 0) {
    join_cond.wait(join_lock);
  }
  join_lock.unlock();
}

size_t fiber_control::get_tid() {
  fiber_control::tls* tls = get_tls_ptr();
  if (tls != NULL) return reinterpret_cast<size_t>(tls->cur_fiber);
  else return (size_t)(0);
}


bool fiber_control::in_fiber() {
  return get_tls_ptr() != NULL;
}

void fiber_control::deschedule_self(pthread_mutex_t* lock) {
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

size_t fiber_control::get_worker_id() {
  fiber_control::tls* tls = get_tls_ptr();
  if (tls != NULL) return tls->workerid;
  else return (size_t)(-1);
}

void fiber_control::schedule_tid(size_t tid, bool priority) {
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
    fib->priority = priority;
    fib->lock.unlock();
    //size_t choice = (fib->affinity >= 0) ? fib->affinity : fib->parent->load_balanced_worker_choice(fib->id);
    size_t choice = 0;
    if (fib->affinity >= 0) choice = fib->affinity;
    else choice = get_worker_id();

    if (choice == (size_t)(-1)) {
      choice = fib->parent->load_balanced_worker_choice(fib->id);
    }
    fib->parent->reschedule_fiber(choice, fib);
  } else {
    //printf("Scheduling requested of running thread %ld\n", fib->id);
    fib->lock.unlock();
  }
}


void fiber_control::set_tls_deleter(void (*deleter)(void*)) {
  flsdeleter = deleter;
}

void* fiber_control::get_tls() {
  fiber_control::tls* f = get_tls_ptr();
  if (f != NULL) {
    return f->cur_fiber->fls;
  } else {
    // cannot get TLS of a non-fiber
    ASSERT_MSG(false, "Trying to get a fiber TLS from a non-fiber");
    return NULL;
  }
}

void fiber_control::set_tls(void* tls) {
  fiber_control::tls* f = get_tls_ptr();
  if (f != NULL) {
    f->cur_fiber->fls = tls;
  } else {
    // cannot get TLS of a non-fiber
    ASSERT_MSG(false, "Trying to get a fiber TLS from a non-fiber");
  }
}


void fiber_control::instance_set_parameters(size_t nworkers = 0,
                                            size_t affinity_base = 0) {
  instance_construct_params_nworkers = nworkers;
  instance_construct_params_affinity_base = affinity_base;
}

fiber_control& fiber_control::get_instance() {
  fiber_control::instance_created = true;
  // set sane defaults
  if (instance_construct_params_nworkers == 0) {
    instance_construct_params_nworkers = thread::cpu_count();
  }
  static fiber_control singleton(instance_construct_params_nworkers, 
                                 instance_construct_params_affinity_base);
  return singleton;
}

}
