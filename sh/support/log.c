/**
 * support/log.[ch]
 * logging, an alternative to d{0,1}printf
 *
 * (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
 * email at willem _-at-_ computer DOT org
 *
 * 3-clause BSD applies
 * */

#ifdef __KERNEL__
#include <linux/kernel.h>
#include <linux/module.h>
#include <linux/string.h>
#else
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "../wrap/file.h"
#endif

#include "log.h"

static char whitespace[] = "            ";
static int sl_loglevel = LOG_STOP;

#ifdef __KERNEL__
#define __print(fd, text) printk("%s", text) 	// weird format because of GCC 4.2 check
#else
#define __print(fd, text) __orig_write(fd, text, strlen(text))
#endif

static void
__write_whitespace(int fd, int len)
{
#ifndef NDEBUG
	if (len > 8) { // hardcoded to be below true length
		__print(fd, "out of whitespace\n");
	}
#endif
	if (len) {
		whitespace[len] = '\0';
		__print(fd, whitespace);
		whitespace[len] = ' ';
	}
}

static void
__write(int fd, enum logtype level, const char *pre, const char *line)
{


	if (level <= sl_loglevel) {
		// write generic header, identifying level, type, etc.
#if !defined __KERNEL__ && !defined NDEBUG 
		{
			char pidbuf[16];
			snprintf(pidbuf, 15, "[%u]", getpid());
			__print(fd, pidbuf);
		}
#endif
		__print(fd, pre);
		
		// add whitespace padding depending on level
		__write_whitespace(fd, level);
		
		// write actual message
		__print(fd, line);
		__print(fd, "\n");
	}
}

#ifdef __KERNEL__

// TODO: move to using our own logging buffer
static void 
write_log(enum logtype level, const char *pre, const char *line)
{
	__write(0, level, pre, line);
}

int 
log_init(void *unused)
{
	return 0;
}

int 
log_exit(void *unused)
{
	return 0;
}

#else

static int logfd = -1;

/** write a message to the log */
static void 
write_log(enum logtype level, const char *pre, const char *line)
{
	if (logfd >= 0) {
		__write(logfd, level, pre, line);
	
		// when debugging, copy important messages to screen
#if !defined NDEBUG 
		if (logfd > 2 && level <= LOG_WARN) 
			__write(2, level, pre, line);
#endif
	}
}

#define MAXNAME 64
static char name[MAXNAME + 1];

int 
log_init(void* unused)
{
	char *tmpdir, *user;
	char linkname[MAXNAME + 1];
	int loglevel_set = 0;

	if (getenv("LOGLEVEL")) {
		sl_loglevel = strtol(getenv("LOGLEVEL"), NULL, 10);
		loglevel_set = 1;
	}

	// log to terminal?
	if (getenv("LOGTERM")) {
		logfd = 1;
		sl_log(LOG_LOW, "logging to terminal");
		return 0;
	}

	// get some metadata to name the file descriptively
	// NB: this is unsafe. check that it is truly a dir?
	tmpdir = getenv("TMPDIR");
	if (!tmpdir)
		tmpdir = "/tmp";
	user = getlogin();
	if (!user)
		user = getenv("USER");

	// create and open a new log file
	snprintf(name, MAXNAME, "%s/streamline.%s.%lu.log", 
		tmpdir, user, time(NULL));
	logfd = __orig_open(name, O_WRONLY | O_CREAT, 0644);
	if (logfd < 0) {
		fprintf(stderr, "error opening log\n");
		return -1;
	}
	
	// set the 'latest' symlink to this file
	snprintf(linkname, MAXNAME, "%s/streamline.%s.latest.log",
		 tmpdir, user);
	unlink(linkname); // don't care whether there was a link before
	if (link(name, linkname) < 0)
		sl_log(LOG_WARN, "error linking log\n");

	// we set this before, but defer output until initialized
	if (loglevel_set)
		sl_log(LOG_LOW, "set loglevel to %d", sl_loglevel);

	return 0;
}

// choose the output file descriptor
void 
log_set(int fd)
{
	logfd = fd;
}

// get the output file descriptor
int 
log_get(void)
{
	return logfd;
}

// get the output filename (if any)
const char * 
log_getname(void)
{
	if (logfd > 2)
		return name;
	else
		return NULL;
}

int 
log_exit(void* unused)
{
	if (logfd >= 0)
		__orig_close(logfd);
	return 0;
}

#endif

void 
sl_log(enum logtype type, const char *format, ...)
{
#define SLLOGSZ 256
	char buf[SLLOGSZ];
	const char *pre;
	va_list ap;

	va_start(ap, format);
	vsnprintf(buf, SLLOGSZ - 1, format, ap);
	va_end(ap);	
	
	switch (type) {
		case LOG_BUG :   pre = "[BUG  ] "; break;
		case LOG_ERR :   pre = "[ERR  ] "; break;
		case LOG_WARN :  pre = "[WARN ] "; break;
		case LOG_MSG :   pre = "[Info ] "; break;
		case LOG_START : pre = "[Start] "; break;
		case LOG_STOP :  pre = "[Stop ] "; break;
		case LOG_LOW :   pre = "[Info ] "; break;
		case LOG_ALL :   pre = "[Info ] "; break;
		default :        pre = "[LOGBUG]"; break;
	}
	write_log(type, pre, buf);
}

#ifdef __KERNEL__
EXPORT_SYMBOL(sl_log);
EXPORT_SYMBOL(log_init);
EXPORT_SYMBOL(log_exit);
#endif

