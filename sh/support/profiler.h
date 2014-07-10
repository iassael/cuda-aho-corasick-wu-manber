/*
	Fairly Fast Packet Filter

	`stealth' profiling headerfile

	Licensed under the GPL
	Copyright (c) Herbert Bos, 2003-2004

	this version was created by Willem de Bruijn (wdebruij_AT_dds.nl), 2004

	\file
	this `class' implements a clockcycle profiler. When the PROFILE
	macro is set, the code will calculate processor cycle counts
	for abitrary program flows by converting the profiler(x),
	profiler_begin(x) and profiler_end(x) macro's into full functions.

	Anytime one of the profiler functions is encountered, data is
	collected and stored for later calculation.

	Note that by not setting the PROFILE macro, the macro's are not
	expanded and therefore the profiler will have no impact on the
	executables' performance.

	[USAGE] define the variable PROFILER_PROC_NAME somewhere to the
	name of the file under /proc that you want to create. If left
	undefined it will default to "ffpf"

	since profiler uses integers to discriminate among classes, I
	suggest you add macro's that expand to unique classkeys either in
	your own code or (if your code will be bundled with this package)
	below (near PROFILER_HOOK) .

	[NB] the kernel's print function, printk, cannot output floating
	point values. Therefore we have resorted to printing the floats
	in another notation, namely as hex integers.  Use this output
	by converting it to the right representation in userspace. For
	this you could use a perl shellscript or something. Currently,
	no such scripts has been written.
*/

#ifndef PROFILE_H
#define PROFILE_H

#ifndef PROFILER_PROC_NAME
#define PROFILER_PROC_NAME "ffpf"
#endif

#define PROFLEN 101	///< number of samples per class
#define PROFWIDTH 9	///< number of classes

/**
	defines to keep track of profile code
	
	you are advised to use these (add yours),
	so that you don't end up with duplicate calls.
*/
#define PROFILER_HOOK		1
#define PROFILER_FILTER		2
#define PROFILER_COPY		3
#define PROFILER_COPY_PKT	4
#define PROFILER_TEST		5
#define PROFILER_BPF_CHECK	6
#define PROFILER_SIGMIN		7
#define PROFILER_SIGMAX		8


/** the profiler routine stores the processor counter
    @param int class. separates streams of statistics.
*/
inline void __internal_profiler(int class);

/** a more explicit version of exec_profiler(..).
	use this function and exec_profiler_end to be
	sure when data collection starts and finishes.
	Consecutive calls to exec_profiler_begin will
	reset the temporary databuffer. The endresult
	is an offset, instead of the raw values. This
	is probably what you want.

	Note that we have no safety checks in place
	for buffer overflows. That's your resposibility.
*/
inline void __internal_profiler_begin(int class);

/** see exec_profiler_begin .*/
inline void __internal_profiler_end(int class);

/** output profiler data. */
inline void __internal_profiler_show(void);

#ifdef __KERNEL__
/** register to procfs. Automatically calls init_profiler (just in case you forget) */
void __internal_profiler_procfs_open(void);
/** unregister from procfs */
void __internal_profiler_procfs_close(void);
#endif /* __KERNEL__ */

#ifdef PROFILE

#define profiler(x) __internal_profiler(x)
#define profiler_begin(x) __internal_profiler_begin(x)
#define profiler_end(x) __internal_profiler_end(x)

#ifdef __KERNEL__
#define profiler_procfs_open() __internal_profiler_procfs_open()
#define profiler_procfs_close() __internal_profiler_procfs_close()
#else
#define profiler_procfs_open()
#define profiler_procfs_close()
#endif /* __KERNEL__ */

#define profiler_init() __internal_profiler_init()
#define profiler_show() __internal_profiler_show()

#else /* PROFILE */

#define profiler(x)
#define profiler_begin(x)
#define profiler_end(x)

#define profiler_procfs_open()
#define profiler_procfs_close()

#define profiler_init()
#define profiler_show()

#endif /* PROFILE */

#endif /* PROFILE_H */

