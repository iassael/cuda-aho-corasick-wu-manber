/**
 * support/filedes.[ch]
 * support incoming signals (such as SIGIO in POSIX userspace)
 *
 * (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
 * email at willem _-at-_ computer DOT org
 *
 * 3-clause BSD applies
 * 
 * Modified by Tudor Zaharia on Aug. 17 2010
 * tudor _at_ microcontroller DOT ro
 * - sinchronized access to slrun on each fd
 * */

#ifdef __KERNEL__
#else
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/time.h>
#include <errno.h>
#include <assert.h>
#endif

#include <semaphore.h>

#include "../support/macros.h"
#include "../support/log.h"
#include "../support/timer.h"
#include "../wrap/file.h"
#include "../wrap/origsocket.h"
#include "../core/datapath.h"
#include "filedes.h"

struct sighandler {
	enum slsig_action action;
	union {
		struct instance *instance;
		void (*callback)(int fd);
	} ptr;
	int backtrack_fd;	/// used to close clientfds after an acceptfd
};

#define MAXFD 512
static struct sighandler * handlerlist[MAXFD];

// these semaphores are used for synchronizing access to the slrun_slow()
static sem_t sems[MAXFD];

/******** implementation-specific support code **********/

int fd_setasync(int fd)
{
#if linux
	int flags;
	
	flags = __orig_fcntl(fd, F_GETFL, 0);
	if (flags < 0)
		goto err;
	if (__orig_fcntl(fd, F_SETFL, flags | O_NONBLOCK | O_ASYNC) == -1)
		goto err;
	if (__orig_fcntl(fd, F_SETOWN, getpid()) == -1)
		goto err;
	return 0;
err:
	sl_log(LOG_WARN, "failed to set fd to async");
	return -1;
#else
	sl_log(LOG_WARN, "async IO not supported"); 
	return -1;
#endif
}

/******** core functions: handlers and callback **********/

/** read out a filedescriptor 
 *
 * because most callbacks will read from a buffer we implemented this
 * functionality locally */
static void sigaction_read(unsigned long sigid, struct instance *instance) 
{
#define MAX_LINESZ 1500 
	char data[MAX_LINESZ];
	int size, total = 0;
	int fd = (int) sigid;

        // enter critical section
        if ( -1 == sem_wait(&sems[fd]) ) perror("semop error");

	size = __orig_read(fd, data, MAX_LINESZ);
	while (size > 0) {
		total += size;
		slrun_slow(instance, data, size);
		size = __orig_read(fd, data, MAX_LINESZ);
	};

	if (size < 0 && errno != EAGAIN)
		perror("read()");
	else if (size == 0) {	// EOF
		slrun_slow(instance, NULL, 0);
	}

	// leave critical section
        if ( -1 == sem_post(&sems[fd]) ) perror("semop error");
}

static void sigaction_accept(unsigned long sigid, struct instance *instance)
{
	int fd = (int) sigid;
	int client_fd;

	client_fd = __orig_accept(fd, NULL, NULL);
	if (client_fd < 0) {
		perror("accept()");
		return;
	}
	filedes_add(client_fd, instance, SIGH_READ);
	handlerlist[client_fd]->backtrack_fd = sigid;
}

/** call a process2() member */
static void sigaction_process(unsigned long sigid, struct instance *instance) 
{
	instance->fdata.func->process2(NULL, NULL, &instance->fdata);
}

/* Handle a SIGIO signal.
 *
 * On receiving one of these two signals, this function
 * executes a non-blocking select() over all file descriptors 
 * registered to support/filedes. On return, it executes all
 * registered handlers for the descriptors on which data is
 * available.
 * */
static void signal_callback(int signal)
{
	struct timeval tv = { .tv_sec = 0, .tv_usec = 0};
	fd_set readfds;
	int i, highest_fd = -1, total;

	// only handle registered signals
	if (signal != SIGIO)
		return;
	
	// add all file descriptors to the listen set
	FD_ZERO(&readfds);
	for(i = 0; i < MAXFD; i++) {
		if (handlerlist[i]) {
			FD_SET(i, &readfds);
			highest_fd = i;
		}
	}
	
	// listen on the descriptor set
	total = __orig_select(++highest_fd, &readfds, NULL, NULL, &tv);
	if (total < 0) {
		if (errno != EINVAL)
			dprintf("error in filedes select\n");
		return;
	}

	// trigger actions for all descriptors on which data is waiting
	for(i = 0; total && i <= highest_fd; i++) {
		if (FD_ISSET(i,&readfds)) {
			switch (handlerlist[i]->action) {
				case SIGH_PROCESS : 
				sigaction_process(i, handlerlist[i]->ptr.instance);
				break;
				case SIGH_READ :
				sigaction_read(i, handlerlist[i]->ptr.instance);
				break;
				case SIGH_ACCEPT :
				sigaction_accept(i, handlerlist[i]->ptr.instance);
				break;
				case SIGH_CALLBACK : 
				handlerlist[i]->ptr.callback(i);
				break;
			};
			total--;
		}
	}
}

/******** bookkeeping **********/

// install the SIGIO handler
int filedes_init(void* unused)
{
	signal(SIGIO, signal_callback);
	return 0;
}

int filedes_exit(void *unused)
{
	signal(SIGIO, SIG_DFL);
	return 0;
}

int filedes_add(int fd, void *ptr, enum slsig_action action)
{
	struct sighandler * sigh;

	assert(fd < MAXFD);
	if (handlerlist[fd])
		return -1;

        // create the semaphore
        if (sem_init(&sems[fd], 0, 1) == -1)
                return -1;
	
	// create the structure
	sigh = myalloc(sizeof(struct sighandler));
	sigh->ptr.instance = (struct instance *) ptr;
	sigh->backtrack_fd = -1;
	sigh->action = action;
	
	// add it to the list
	handlerlist[fd] = sigh;

	// ask the OS to signal us when data arrives on this fd. 
	fd_setasync(fd);
	
	// bootstrap first read (for files)
	if (action == SIGH_READ || action == SIGH_PROCESS)
		signal_callback(SIGIO);

	return 0;
}

/** close all connections that depend on the parameter */
static void __signalhandler_deldep(int fd)
{
	int i;

	for(i=0; i<MAXFD; i++)
		if (handlerlist[i] &&
		    handlerlist[i]->backtrack_fd == fd) {
			__orig_close(i);
			myfree(handlerlist[i]);
			handlerlist[i] = NULL;
		}
}

/** remove a connection. may also recursively remove dependent connections */
int filedes_del(int fd)
{
	struct sighandler * sigh;
	
	assert(fd < MAXFD);
	if (!handlerlist[fd])
		return -1;

	sigh = handlerlist[fd];
	handlerlist[fd] = NULL;

	// close all client connections if this is an accept descriptor.
	if (sigh->action == SIGH_ACCEPT)
		__signalhandler_deldep(fd);
	
	myfree(sigh);

        // destroy semaphore
        sem_destroy(&sems[fd]);

	return 0;
}

