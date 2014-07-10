// lock.h
// mutual exclusion and other locking support
//
// (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// LGPL license applies

// this is a wrapper around various locking methods
// note: slmutex_trylock returns !0 if a lock is held, 0 otherwise


#ifdef __KERNEL__
#include <linux/version.h>

#if LINUX_VERSION_CODE < KERNEL_VERSION(2,6,17)

#include <asm/semaphore.h>
#define slmutex struct semaphore
#define slmutex_init(my_mutex) init_MUTEX(my_mutex)
#define slmutex_static(my_mutex) DECLARE_MUTEX(my_mutex)
#define slmutex_extern(my_mutex) extern struct semaphore my_mutex
#define slmutex_lock(my_mutex) do {} while (down_interruptible(my_mutex))
#define slmutex_unlock up
#define slmutex_trylock(my_mutex) (down_trylock(my_mutex) ? 0 : 1)

#else // newer kernel?

#include <linux/mutex.h>
#include <linux/hardirq.h>
#include "macros.h"
#define slmutex struct mutex
#define slmutex_init(my_mutex) mutex_init(my_mutex)
#define slmutex_static DEFINE_MUTEX
#define slmutex_extern(my_mutex) extern struct mutex my_mutex
#define slmutex_trylock mutex_trylock
#if 1
#define slmutex_lock mutex_lock
#define slmutex_unlock mutex_unlock
#else
#define slmutex_lock(my_mutex) 						      \
	do {dprintf("mutex_lock in %s. atomic=%c locked=%c\n", __FUNCTION__,  \
		    in_atomic()?'y':'n', mutex_is_locked(my_mutex)?'y':'n');  \
	    mutex_lock(my_mutex); 					      \
	    dprintf("mutex locked\n");					      \
	} while(0)
#define slmutex_unlock(my_mutex) 					      \
	do {dprintf("mutex_unlock in %s. atomic=%c locked=%c\n", __FUNCTION__,\
		    in_atomic()?'y':'n', mutex_is_locked(my_mutex)?'y':'n');  \
	    mutex_unlock(my_mutex); 					      \
	    dprintf("mutex_unlocked\n"); 				      \
	} while(0)
#endif
#endif

#else // userspace?

#define in_atomic() (0)

#include <pthread.h>
#define slmutex pthread_mutex_t
#define slmutex_init(my_mutex) pthread_mutex_init(my_mutex, NULL);
#define slmutex_static(my_mutex) pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER
#define slmutex_extern(my_mutex) extern pthread_mutex_t my_mutex
#define slmutex_trylock(my_mutex) (pthread_mutex_trylock(my_mutex) ? 0 : 1)
#if 1
#define slmutex_lock pthread_mutex_lock
#define slmutex_unlock pthread_mutex_unlock
#else
#define slmutex_lock(my_mutex) 						      \
	do {dprintf("mutex_lock %p in %s\n", my_mutex, __FUNCTION__);         \
	    pthread_mutex_lock(my_mutex); 		                      \
	    dprintf("mutex locked\n");					      \
	} while(0)
#define slmutex_unlock(my_mutex) 					      \
	do {dprintf("mutex_unlock %p in %s\n", my_mutex, __FUNCTION__);       \
	    pthread_mutex_unlock(my_mutex); 				      \
	} while(0)
#endif
#endif

