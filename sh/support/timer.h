// timer.[ch]
// wrapper around OS-specific alarm signals
//
// (c) 2005, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_AT_- computer.org
//
// LGPL License applies

#ifndef SL_SUPPORT_TIMER
#define SL_SUPPORT_TIMER

#include <streamline/signal.h>

/** execute a task in the background
 *
 *  @param recur sets how often the task should be executed, 
 *         -1 for indefinite or until task_stop is called.
 *
 *  @return an opaque pointer to pass to task_stop */
void * task_start(void(*func)(void*), void * arg, long recur, long timeout);
void task_stop(void *task);

/** Cancel all outstanding tasks.
 *  Some tasks may still fire, but all are stopped
 *  when this function returns.
 *
 *  return 0 on success, failure otherwise */
int task_stop_all(void);


#endif

