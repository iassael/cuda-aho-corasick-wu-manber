// timer.[ch]
// wrapper around OS-specific alarm signals
//
// (c) 2005, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_AT_- computer.org
//
// LGPL License applies

#ifdef __KERNEL__
#include <linux/kernel.h>
#include <linux/module.h>
#include <linux/version.h>
#include <linux/wait.h>
#include <linux/sched.h>
#include <linux/workqueue.h>
#include <linux/limits.h>
#include <asm/atomic.h>
#else
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <signal.h>
#include <unistd.h>
#include <setjmp.h>
#include <pthread.h>
#include <limits.h>
#include "../wrap/file.h"
#endif
#include <linux/param.h>

#include "../core/config.h"
#include "list.h"
#include "log.h"
#include "macros.h"
#include "timer.h"
#include "lock.h"

#ifdef __KERNEL__

/// the number of active tasks. 
//  is forced to 0 by interrupt_deep to cancel all tasks
int tasks_stop = 0;

struct list * timers;

struct task {
	unsigned long jiffies;
	unsigned long recur;
	void (*func)(void *);
	void * arg;
	int forced_stop;
#if LINUX_VERSION_CODE < KERNEL_VERSION(2,6,20)
	struct work_struct dws;
#else
	struct delayed_work dws;
#endif
};

// callback: calls the function and reenables the timer
#if LINUX_VERSION_CODE < KERNEL_VERSION(2,6,20)
static void task_callback(void * arg) {
	struct task * t = arg;
#else
static void task_callback(struct work_struct * ws) {
	struct delayed_work * dw = container_of(ws, struct delayed_work, work);
	struct task * t = container_of(dw, struct task, dws);
#endif
	
	if (!tasks_stop && !t->forced_stop) {
		if (t->recur > 0)
			t->recur--;
		
		if (unlikely(!t->func)) 
			dprintf("ERR at %s.%d", __FUNCTION__, __LINE__);
		else	
			t->func(t->arg);
		
		if (t->recur) {
			schedule_delayed_work(&t->dws, t->jiffies);
			return;
		}
	}

	kfree(t);
}

void * task_start(void(*func)(void*), void * arg, long recur, long timeout)
{
	struct task * t;

	// fill our structure
	t = kzalloc(sizeof(struct task), GFP_ATOMIC);
	if (!t) {
		sl_log(LOG_ERR, "out of atomic memory");
		return NULL;
	}
	t->func = func;
	t->arg = arg;
	t->recur = recur;
	t->jiffies = (HZ * timeout) / 1000000;

	// initialize the waitqueue element
#if LINUX_VERSION_CODE < KERNEL_VERSION(2,6,20)
	INIT_WORK(&t->dws, task_callback, t);
#else
	INIT_DELAYED_WORK(&t->dws, task_callback);
#endif
	
	if (t->jiffies)
		schedule_delayed_work(&t->dws, t->jiffies);
	else
#if LINUX_VERSION_CODE < KERNEL_VERSION(2,6,20)
		schedule_work(&t->dws);
#else
		schedule_work(&t->dws.work);
#endif
	return t;
}

void task_stop(void *task)
{
	struct task * t = task;

	t->forced_stop = 1;
	cancel_delayed_work(&t->dws);
	flush_scheduled_work();
}

int task_stop_all(void)
{
	tasks_stop = 1;
	flush_scheduled_work();
	return 0;
}

int usleep_deep(long timeout)
{
	set_current_state(TASK_INTERRUPTIBLE);
	schedule_timeout(usecs_to_jiffies(timeout));
	return 0;
}

int interrupt_deep(void)
{
	dprintf("%s called in kernel: unsupported\n", __FUNCTION__);
	return -1;
}

EXPORT_SYMBOL(task_start);
EXPORT_SYMBOL(task_stop);
EXPORT_SYMBOL(task_stop_all);

EXPORT_SYMBOL(interrupt_deep);
EXPORT_SYMBOL(usleep_deep);

#else /* !__KERNEL__ */

struct task_args {
	void (*func)(void*);
	pthread_t thread;
	void * arg;
	long timeout;
	long recur;
};

// TODO: use atomic types
static int tasks_stop;
static int tasks_active; 

static void * 
delayed_func(void *thread_arg) 
{
	struct task_args *ta = thread_arg;

	while (ta->recur > 0 || ta->recur == -1) {
		if (usleep_deep(ta->timeout) < 0) {
			ta->recur = 0;
			break;
		}
		
		if (tasks_stop)
			break;

		ta->func(ta->arg);

		if (ta->recur > 0)
			ta->recur--;
	}
	
	// task_stop_all does not call pthread_join
	// and noone is waiting if recurrence ended
	if (tasks_stop || ta->recur == 0) {
		pthread_detach(ta->thread);
		myfree(ta);
		tasks_active--;
	}

	return NULL;
}

void * task_start(void(*func)(void*), void * arg, long recur, long timeout)
{
	struct task_args *ta;

	if (tasks_stop)
		return NULL;

	ta = mycalloc(1, sizeof(struct task_args));
	ta->func = func;
	ta->arg = arg;
	ta->timeout = timeout;
	ta->recur = recur;
	tasks_active++;

	pthread_create(&ta->thread, NULL, delayed_func, ta);
	return ta;
}

// don't allow the purging of all tasks interfere with a single task
// that is to be closed
slmutex_static(task_mutex);

void task_stop(void *task)
{
	struct task_args *ta = task;

	slmutex_lock(&task_mutex);
	if (ta) {
		if (!tasks_active)
			sl_log(LOG_BUG, "waiting for nonexistent task");
		ta->recur = -2; // stop, signal that we will wait for the result
		interrupt_deep();
		pthread_join(ta->thread, NULL); 
		myfree(ta);
		tasks_active--;
	}
	slmutex_unlock(&task_mutex);
}

/// May only be called from process context, because it may sleep.
int task_stop_all(void)
{
	slmutex_lock(&task_mutex);
	if (tasks_active) {
		tasks_stop = 1;
		if (interrupt_deep()) {
			sl_log(LOG_ERR, "Failed to interrupt");
			tasks_stop = 0;
			tasks_active = 0; // try to set to a 'stable' state
			return -1;
		}
		while (tasks_active) {
			sl_log(LOG_MSG, "Waiting for %d tasks to finish\n", tasks_active);
			sleep(1);
			tasks_stop = 0;
		}
	}
	slmutex_unlock(&task_mutex);
	
	return 0;
}

// HACKHACKHACK replace with nice open on load + close on unload
static int shallowfd = -1;
static int deepfd = -1;

static int
__usleep_sl(int *fd, const char *name, long timeout)
{
	if (unlikely((*fd) == -1)) {
		(*fd) = __orig_open(name, O_WRONLY);
		if ((*fd) < 0) {
			sl_log(LOG_ERR, "open timer failure");
			return -1;
		}
	}
	return __orig_write((*fd), &timeout, sizeof(long));
}

// pause the thread for the given number of microseconds
//
// We try to avoid having to use POSIX signals. If kernelspace Streamline
// exists, we use the sysfs timer file, otherwise we rever to SIGALRM
//
// returns 0 on success. timeout left if > 0, signal arrived if < 0
int usleep_deep(long timeout)
{
	return __usleep_sl(&deepfd, SYSFS_TIMER_DEEP, timeout);
}

// pause the thread for the given number of microseconds
// or until a streamline signal arrives
//
// returns 0 on success. timeout left if > 0, signal arrived if < 0
int usleep_shallow(long timeout)
{
	return __usleep_sl(&shallowfd, SYSFS_TIMER_SHALLOW, timeout);
}

int interrupt_deep(void)
{
	int fd;
	char useless = 0;

	fd = __orig_open(SYSFS_TIMER_INTERRUPT, O_WRONLY);
	if (unlikely(fd < 0)) {
		sl_log(LOG_LOW, "failed to call deep interrupt. POSIX timers?");
		return -1;
	}
	
	if (unlikely(__orig_write(fd, &useless, 1)) < 0)
		return -1;

	if (unlikely(__orig_close(fd)))
		return -1;
	
	return 0;
}

#endif /* !__KERNEL__ */

