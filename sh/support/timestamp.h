// timestamp.h
// location independent timestamping
//
// (c) 2005, willem de bruijn, vrije universiteit amsterdam
// email at willem -_at_- computer.org
//
// BSD license applies


#ifdef __KERNEL__
#include <linux/kernel.h>
#include <linux/version.h>
#if LINUX_VERSION_CODE < KERNEL_VERSION(2,6,18)
#include <linux/config.h>		// TSC available?	
#endif
#include <linux/time.h>
#include <linux/timex.h>		// platform independent backup
#ifdef CONFIG_X86_TSC
#include <linux/cpufreq.h>		// cpufreq. a lousy method	
#include <asm/msr.h>			// 64bit cycle-accurate counter
#include <linux/types.h>
#endif 
#else
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#endif

#ifdef CONFIG_ARM
/// hack TODO: fix
#define cpu_khz 600000
#endif

// timestamp_get
#if (defined __KERNEL__ && defined CONFIG_X86_TSC) || !defined NO_X86
typedef uint64_t tstamp_t;
static inline uint64_t timestamp_get(void) {
	register uint32_t count_low, count_high;
	asm("rdtsc" :"=a" (count_low), "=d" (count_high));
	return ( ((uint64_t) count_high) << 32) + count_low;
}
#else
#ifdef __KERNEL__
typedef cycles_t tstamp_t;
#define timestamp_get get_cycles
#else
typedef clock_t tstamp_t;
#define timestamp_get clock
#endif
#endif

// timestamp_to
#ifdef __KERNEL__
static inline tstamp_t timestamp_to(int sec, int usec)
{
	return (cpu_khz * usec) + (cpu_khz * 1000 * sec);
}
#else
static inline tstamp_t timestamp_to(int sec, int usec)
{
	return (CLOCKS_PER_SEC * sec) + ((CLOCKS_PER_SEC/1000000) * usec);
}
#endif 


