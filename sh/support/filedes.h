/**
 * support/filedes.[ch]
 * support incoming signals (such as SIGIO in POSIX userspace)
 *
 * (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
 * email at willem _-at-_ computer DOT org
 *
 * 3-clause BSD applies
 * */

#ifndef SL_SUPPORT_FILEDES_H
#define SL_SUPPORT_FILEDES_H

enum slsig_action {SIGH_PROCESS=1, SIGH_READ, SIGH_ACCEPT, SIGH_CALLBACK};

int filedes_add(int fd, void *ptr, enum slsig_action action);
int filedes_del(int fd);

int filedes_init(void*);
int filedes_exit(void*);

#endif

