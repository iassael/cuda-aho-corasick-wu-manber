// prettyprint.[ch]
// print non-trivial data to screen
//
// (c) 2007, Willem de Bruijn, Vrije Universiteit Amsterdam
// email at willem -_at_- computer.org
//
// LGPL license applies

#ifdef __KERNEL__
#include <linux/kernel.h>
#include <linux/module.h>
#include <linux/if_ether.h>	
#include <linux/in.h>	
#include <linux/ip.h>
#include <linux/tcp.h>
#include <linux/udp.h>
#include <linux/icmp.h>
#else
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#if linux
#include <netinet/if_ether.h>
#include <netinet/ip.h>
#include <netinet/tcp.h>
#include <netinet/udp.h>
#include <netinet/ip_icmp.h>
#else
#include "../hw/proto.h"
#endif
#endif

#include "macros.h"
#include "lock.h"
#include "prettyprint.h"

/** Pretty print data in vertically split [hex | decimal] notation 
 *
 *  The passed string must be at least 80 bytes.
 * */
int
writedata(char *out, int olen, const char *data, int dlen) {
    	int i = 0, off = 0;
	int elem;

	olen--; // leave room for the terminating '\0'
	do {
		// phase 1: print hex
		for (elem = 0; elem < HEXWIDTH && i + elem < dlen; elem++) {
			off += snprintf(out + off, olen - off, "%x%x ", 
				        (data[i + elem] & 0xf0) >> 4,
				         data[i + elem] & 0xf);
			
			if (elem == (HEXWIDTH / 2) - 1) {
				out[off] = ' ';
				off++;
			}
		}
		
		// fill out the last line
		for (; elem < HEXWIDTH; elem ++) {
			out[off] = ' ';
			out[off + 1] = ' ';
			out[off + 2] = ' ';
			off += 3;
		}
		
		// insert room
		out[off] = ' ';
		out[off + 1] = ' ';
		out[off + 2] = ' ';
		off += 3;

		// phase 2: print visible
		for (elem = 0; elem < HEXWIDTH && i + elem < dlen; elem++) {
			if (data[i + elem] >= 32 && data[i + elem] < 126)
				out[off + elem] = data[i + elem];
			else
				out[off + elem] = '.';
		}
		off += elem;
		out[off] = '\n';
		off++;
		i += HEXWIDTH;
	} while (i < dlen && off < olen);

	out[off] = '\n';
	off++;
	out[off] = '\0';
	return off;
}

void
displaydata(const char *data, int dlen)
{
	char *out;
	int len, mlen;
#ifndef __KERNEL__
	int ret;
#endif
	
	if (dlen) {
		len = 5 * dlen;
		mlen = max(len + 1, 800);

		// allocate the block. difficult only because of 
		// possible execution in kernel interrupt context.
#ifdef __KERNEL__
		if (in_interrupt()) {
			out = kmalloc(mlen, GFP_ATOMIC);
			memset(out, 0, mlen);
		}
		else
#endif		
			out = mycalloc(mlen, 1);

		// malloc failed error handling
		if (!out) {
			const char error[] = "BUG in displaydata\n";
#ifdef __KERNEL__
			printk(error);
#else
			ret = write(1, error, 20);
#endif
			return;
		}

		// fill in contents and write output
		len = writedata(out, len - 1, data, dlen);
		out[len] = '\0';
#ifdef __KERNEL__
		printk("%s", out);
#else
		ret = write(1, out, len + 1);
#endif
		myfree(out);
	}
}

/** Prettyprint an IP address.
 *  @returns the number of bytes written */
int
writeip(char * data, int dlen, const uint8_t* ip, uint16_t port)
{
	int res; 
#ifdef __KERNEL__
	res = snprintf(data, dlen, "%hu.%hu.%hu.%hu", ip[0], ip[1], ip[2], ip[3]);
#else
	res = snprintf(data, dlen, "%hhu.%hhu.%hhu.%hhu", ip[0], ip[1], ip[2], ip[3]);
#endif
	if (port)
		res += snprintf(data + res, dlen - res, ":%hu", ntohs(port));
	return res;
}

/** Print an ip address to stdout (w/o endline) */
void 
displayip(const uint8_t* ip, uint16_t port)
{
	char buf[25];
	writeip(buf, 24, ip, port);
	aprintf("%s", buf);
}

int
writepktinfo(char *out, int olen, const char *pkt, unsigned int plen)
{
	const struct ethhdr *eth = (struct ethhdr *) pkt;
	uint16_t sport=0, dport=0, off, i;
	
	olen--; // leave room for the terminating '\0'
	if (plen < ETH_HLEN)
		return snprintf(out, olen, "%dB: too small for ethernet\n", plen);

	off = snprintf(out, olen, "eth(type %hx, src ", ntohs(eth->h_proto));
	for (i = 0; i < 6; i++)
		off += snprintf(out + off, olen - off, "%hx%hx.", 
				(eth->h_source[i] & 0xf0) >> 4, 
				 eth->h_source[i] & 0xf);
	off += snprintf(out + off, olen - off, ", dest ");
	for (i = 0; i < 6; i++)
		off += snprintf(out + off, olen - off, "%hx%hx.", 
			 	(eth->h_dest[i] & 0xf0) >> 4, 
				 eth->h_dest[i] & 0xf);
	off += snprintf(out + off, olen - off, ")\n");
	
	if ((uint16_t) ntohs(eth->h_proto) == ETH_P_IP) { 
		const struct iphdr *iph;
		
		iph = (struct iphdr*) (pkt + ETH_HLEN);
		off += snprintf(out + off, olen - off, 
				"ip (proto %hu, ttl %hu, ihl %hu, total_len %hu,"
				" src %hu.%hu.%hu.%hu, dst %hu.%hu.%hu.%hu)\n", 
			        iph->protocol,
			        iph->ttl,
			        iph->ihl,
			        ntohs(iph->tot_len),
			        iph->saddr & 0xff, 
			        (iph->saddr & 0xff00) >> 8, 
			        (iph->saddr & 0xff0000) >> 16, 
			        (iph->saddr &0xff000000) >> 24,
			        iph->daddr & 0xff, 
			        (iph->daddr & 0xff00) >> 8, 
			        (iph->daddr & 0xff0000) >> 16, 
			        (iph->daddr &0xff000000) >> 24);
		
		if (iph->protocol == 6 /* TCP */){
			// start of udp and tcp headers are identical. but the following is
			// a bit hackish, I admit
			const struct tcphdr *tcph;
			
			tcph = (struct tcphdr*) ( ((unsigned long) iph) + (4 * ((char) iph->ihl)));
			sport = tcph->source;
			dport = tcph->dest;
			off += snprintf(out + off, olen - off, 
					"tcp(len=%u seqno=%u)\n", 
				        ntohs(iph->tot_len), ntohl(tcph->seq));
		} 
		else if (iph->protocol == 17 /* UDP */){
			const struct udphdr *trans;
			
			trans = (struct udphdr*) ( ((unsigned long) iph) + (4 * ((char) iph->ihl)));
			sport = trans->source;
			dport = trans->dest;
			off += snprintf(out + off, olen - off, "udp(len=%u)\n", 
					ntohs(iph->tot_len));
		}
		else if (iph->protocol == 1 /* ICMP */){
			const struct icmphdr *icmph;
			
			icmph = (struct icmphdr*) ( ((unsigned long) iph) + (4 * ((char) iph->ihl)));
			off += snprintf(out + off, olen - off, 
					"icmp(len=%u type=%hu seq=%hu)\n", 
				        ntohs(iph->tot_len), icmph->type, 
					icmph->un.echo.sequence);
		}
	}
	else
		off += snprintf(out + off, olen - off, "unknown()\n");

	out[off] = '\n';
	off++;
	out[off] = '\0';
	return off;
}

void
displaypktinfo(const void *data, int len)
{
	char *out; 
#ifndef __KERNEL__
	int ret;
#endif

	out = myalloc(240); // 3 lines is the upper limit
	len = writepktinfo(out, 239, data, len);
	out[len] = '\0';
#ifdef __KERNEL__
	printk("%s", out);
#else
	ret = write(1, out, len + 1);
#endif
	myfree(out);
}

