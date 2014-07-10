// macros.[ch]
// simple macros that I reuse often
//
// (c) 2005, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_AT_- computer.org
//
// BSD License applies

#ifndef WJDB_SUPPORT_H
#define WJDB_SUPPORT_H

///////////// KERNELSPACE/USERSPACE COMPAT
#ifdef __KERNEL__
#include <linux/kernel.h>
#include <linux/sched.h>
#include <linux/version.h>
#include <linux/slab.h>
#define myalloc(block) kmalloc(block, GFP_KERNEL)
#define myfree(block) kfree(block)
static inline void * mycalloc(size_t nmemb, size_t size) 
{
	void *data = kmalloc(nmemb * size, GFP_KERNEL);
	if (data)
		memset(data, 0, nmemb * size);
	return data;
}
#define clock() get_cycles()
#define MY_CLOCKRATE 1800000000
#define my_gettimeofday do_gettimeofday
#if LINUX_VERSION_CODE >= KERNEL_VERSION(2, 6, 29)
#define getpid() (current->pid)
#define getuid() (current_uid())
#define getgid() (current_gid())
#else
#define getpid() (current->pid)
#define getuid() (current->uid)
#define getgid() (current->gid)
#endif
#else
#include <stdio.h>
#define myalloc malloc
#define mycalloc calloc
#define myfree(a) free(a)
#define MY_CLOCKRATE CLOCKS_PER_SEC
#define my_gettimeofday(a) gettimeofday(a, NULL)

// PAGE_SIZE is not defined in userspace
#if defined i386 || defined __x86_64__
#define PAGE_SIZE 4096
#else
#define PAGE_SIZE getpagesize()
#endif
// (un)likely is not defined in userspace
#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)
#endif /* __KERNEL__ */

///////////// 32/64 bit COMPAT
#ifdef __LP64__
#define FMT64 "ld"
#define FMT64U "lu"
#else
#define FMT64 "lld"
#define FMT64U "llu"
#endif

///////////// Portable ASSERTIONS
//
// an alternative to assert() that can also work in kernelspace
// NB: it sins against the rule that no control-flow should be in macros...
//     ... but then, so does assert.
// Update (17032008): this stuff is outdated, but 
//     I'm too lazy to clean up all source
#ifdef NDEBUG
#define __check(expression, exec_stmt) \
	do {\
		if (unlikely((expression) == 0)) { \
			exec_stmt; \
		} \
	} while(0)
#else
#define __check(expression,exec_stmt) \
	do { \
		if (unlikely((expression) == 0)) { \
			dprintf("ASSERT FAILED at %s.%d\n",__FUNCTION__,__LINE__); \
			exec_stmt; \
		} \
	} while(0)
#endif

#define check_noop(expression) __check(expression,)
#define check(expression) __check(expression, return -1)
#define check_ptr(expression) __check(expression, return NULL)
#define check_void(expression) __check(expression, return )
#define check_goto(expression) __check(expression, goto cleanup)

#ifdef __KERNEL__
#define assert(stmt) do {if (!(stmt)) panic("ASSERT FAILED at %s.%d\n", __FILE__, __LINE__); } while (0)
#endif

///////////// Portable PRINT
#ifdef __KERNEL__
#define aprintf printk 
#define dprintf printk
#else
#define aprintf printf
#ifdef NDEBUG
#define dprintf(...) 
#else
#define dprintf aprintf
#endif /* NDEBUG */
#endif /* __KERNEL__ */

///////////// Other
#define bug() dprintf("(BUG) at %s:%d\n",__FUNCTION__,__LINE__)
#define returnbug(a) do {bug(); return(a);} while (0)

#ifndef min
#define min(a,b) ( ((a) < (b)) ? (a) : (b) )
#define max(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif

#define is_power2(x) (!(x & (x-1)))

#define __OFF(a, b) (((unsigned long) a) - ((unsigned long) b))

#define IO_IN 1
#define IO_OUT 2

#endif /* WJDB_SUPPORT_H */

