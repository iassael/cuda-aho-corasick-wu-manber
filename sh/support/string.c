// string.c
// standard string functionality that is not always available
//
// (c) 2008, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// LGPL license applies

#ifdef __KERNEL__
#else
#include <stdlib.h>
#endif

#include "macros.h"
#include "log.h"
#include "string.h"

#ifdef __KERNEL__
long 
strtol(const char *in, char **out, int base)
{
	long total=0;
	int i=0;

	if (out || base != 10) {
		sl_log(LOG_ERR, "incomplete strtol called in unsupported mode");
		return 0;
	}

	while (in[i] >= '0' && in[i] <= '9') {
		total *= 10;
		total += in[i] - '0';
		i++;
	}
	return total;
}

// yes, this is an almost exact copy of above. I should've used ##
unsigned long
strtoul_ex(const char *in, char **out, int base, int *err)
{
	unsigned long total=0;
	int i=0;

	if (out || base != 10)
		return 0;

	while (in[i] >= '0' && in[i] <= '9') {
		total *= 10;
		total += in[i] - '0';
		i++;
	}

	// set error if non-digit characters were encountered
	if (err) {
		if (in[i] == '\0')
			*err = 0;
		else
			*err = 1;
	}

	return total;
}

unsigned long 
strtoul(const char *in, char **out, int base)
{
	return strtoul_ex(in, out, base, NULL);
}

char * 
strdup(const char *in)
{
	char *out;
   	int len;

	len = strlen(in);
	out = myalloc(len + 1);
	if (out)
		memcpy(out, in, len);
	out[len]='\0';
	return out;
}
#endif

uint32_t 
strtohost(const char *string, uint16_t *port)
{
	const char * token;
	unsigned short ipseg[4];

	token = strchr(string, ':');
	if (token)
		sscanf(string, "%hu.%hu.%hu.%hu:%hu", 
		       &ipseg[0], &ipseg[1], &ipseg[2], &ipseg[3], port);
	else {
		sscanf(string, "%hu.%hu.%hu.%hu",
		       &ipseg[0], &ipseg[1], &ipseg[2], &ipseg[3]);
		*port = 0;
	}

	return (ipseg[0] << 24) + (ipseg[1] << 16) + (ipseg[2] << 8) 
		+ ipseg[3];
}
