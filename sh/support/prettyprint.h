// prettyprint.[ch]
// print non-trivial data to screen
//
// (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// LGPL license applies

#ifndef __SL_SUPPORT_PRETTYPRINT_H
#define __SL_SUPPORT_PRETTYPRINT_H

#define HEXWIDTH 16

void displaydata(const char *data, int dlen);
void displayip(const uint8_t* ip, uint16_t port);
void displaypktinfo(const void *data, int len);

int writeip(char * data, int dlen, const uint8_t* ip, uint16_t port);
int writedata(char *out, int olen, const char *data, int dlen);
int writepktinfo(char *out, int olen, const char *pkt, unsigned int plen);

#endif

