/**
 * support/log.[ch]
 * logging, an alternative to d{0,1}printf
 *
 * (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
 * email at willem _-at-_ computer DOT org
 *
 * 3-clause BSD applies
 * */

#ifndef SL_SUPPORT_LOG_H
#define SL_SUPPORT_LOG_H

enum logtype {LOG_BUG = 0, LOG_ERR, LOG_WARN, LOG_MSG, LOG_START, LOG_STOP, LOG_LOW, LOG_ALL};

void sl_log(enum logtype type, const char *format, ...);

int log_init(void* unused);
int log_exit(void* unused);

void log_set(int fd);
int log_get(void);
const char * log_getname(void);

#endif

