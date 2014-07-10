/*
	Fairly Fast Packet Filter

	`stealth' profiling sourcefile

	Licensed under the GPL
	Copyright (c) Herbert Bos, 2003-2004

	this version was created by Willem de Bruijn (wdebruij_AT_liacs.nl), 2004
*/

#ifdef PROFILE

#ifdef __KERNEL__
#include <linux/kernel.h>
#include <linux/proc_fs.h>
#else
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CAN_PRINT_FLOATS
#endif

#include "macros.h"
#include "timestamp.h"
#include "profiler.h"

struct profdata {
	tstamp_t cycles[PROFWIDTH][PROFLEN];
	uint32_t index[PROFWIDTH];
};

static struct profdata prof;

/** save a new processor count. */
inline void __internal_profiler(int class){
	/* DIRTY : the first element isn't a diff. 
	 * Either forget about Avg and use only Mean,
	 * or wait long enough for this element to be overwritten. */
	prof.cycles[class][ prof.index[class] ] = timestamp_get();
	prof.index[class] = (prof.index[class] + 1) % PROFLEN;
}

/** start a new processor count calculation. */
inline void __internal_profiler_begin(int class){
	prof.cycles[class][prof.index[class] % PROFLEN] = timestamp_get();
}

/** close a processor count calculation. */
inline void __internal_profiler_end(int class){
	prof.cycles[class][prof.index[class]] = timestamp_get() - prof.cycles[class][prof.index[class]];
	prof.index[class] = (prof.index[class] + 1) % PROFLEN;
	/* note that this would result in a negative result on signed values. We'll have to swap this when calculating results */
}

/* quicksort implementation from wikipedia.org.

   we could have used the qsort(..) function call in userspace, but
   for simplicity we'll use this less optimal algorithm in both kernel-
   and userspace.
*/
void __qsort(tstamp_t* low, tstamp_t* high)
{
   /* We naively use the first value in the array as the pivot */
   /* this will not give good performance real usage */

   tstamp_t * lowbound = low + 1;       /* the high boundary of the low subarray */
   tstamp_t * highbound = high - 1;     /* the low boundary of the high subarray */
   tstamp_t temp;

   while(lowbound <= highbound) /* partition the array */
   {
      if(*lowbound < *low)         /* compare to pivot */
        lowbound++;                /* move lowbound toward the middle */
      else
      {
         temp = *lowbound;         /* swap *lowbound and *highbound */
         *lowbound = *highbound;
         *highbound = temp;
         highbound--;              /* move highbound toward the middle */
      }
   }

   highbound++;                    /* move bounds back to the correct positions */
   lowbound--;

   temp = *low;                    /* move the pivot into the middle */
   *low = *lowbound;
   *lowbound = temp;

   if(low != lowbound)             /* recurse on the subarrays */
     __qsort(low, lowbound);
   if(high != highbound)
     __qsort(highbound, high);
}

tstamp_t __median(int start, int stop, tstamp_t* ldList){
	int middle_floor = start+(stop-start)/2;
	if ( (((stop-start) % 2) + 1) == 1)	// odd number of elements
		return ldList[middle_floor];
	else
		return ((tstamp_t) ( ldList[middle_floor] + ldList[middle_floor + 1]) ) / 2;
}

/**
 * calculate the mean and output information to the standard output queue. 
 * this function is very similar to the one that outputs to procfs.
 * I currently don't have the time to properly merge the two. */
void __internal_profiler_show(void){
	int i, j;
	tstamp_t Q1, Q2, Q3;

	for (i = 0; i < PROFWIDTH; i++) {
		/* are we using this class? */
		if (prof.cycles[i][0]) {
			int val;
			double average;
			
			/* find the last used element in the list */
			val = 0;
			while(val < PROFLEN && prof.cycles[i][val])
				val++;
			if (prof.cycles[i][val])
				val++;

			/* calculate the mean (and the lower and upper quartile) */
			__qsort(&prof.cycles[i][0], &prof.cycles[i][val-1]);
			Q2 = __median(0,val-1, prof.cycles[i]);
			if (val % 2) {
				Q1 = __median(0, val/2 - 1, prof.cycles[i]);
				Q3 = __median(val/2 + 1, val - 1, prof.cycles[i]);
			}
			else{
				Q1 = __median(1, val/2 - 2 , prof.cycles[i]);
				Q3 = __median(val/2 + 1, val - 1, prof.cycles[i]);
			}
			dprintf("class %d: Q1=%llu Q2=%llu Q3=%llu \n", 
				i, Q1, Q2, Q3);

			/* calculate the average */
			average = 0;
			for (j=0; j < val; j++)
				average += ((double) prof.cycles[i][j]) / val;
			dprintf("class %d: average is %lf\n",i, average);
		}
	}
}

#ifdef __KERNEL__ /* we can only export to procfs from the kernel, naturally */

/** static buffer for keeping our fake procfs */
static char procfs_buffer[80 + 5*80 * PROFWIDTH];	/** used for exporting information to procfs */

/** export information to procfs. */
int __internal_profiler_procfs(char *buffer, char **buffer_location, off_t offset, int buffer_length, int zero){
	int len;
	int i, j, val;
	tstamp_t Q1, Q2, Q3;
	double average;

	if (offset > 0)
		return 0;

	memset(procfs_buffer,0, 80 + (5*80 * PROFWIDTH) - 1);
	len = snprintf(procfs_buffer, 17, "kernel profiler\n\n");
	/* Fill the buffer and get its length */
	for (i = 0; i < PROFWIDTH; i++){
		if (prof.cycles[i][0]){	/* are we using this class? then a 0 value is highly unlikely */
			/* find the last used element in the list (might well be PROFLEN */
			val=0;
			while(val < PROFLEN && prof.cycles[i][val]){
				//dprintf("%d,%d:%llu\n",i,val,prof.cycles[i][val]);
				val++;
			}
			if (prof.cycles[i][val])
				val++;

			/* calculate the mean (and the lower and upper quartile) */
			__qsort(&prof.cycles[i][0],&prof.cycles[i][val-1]);
			Q2 = __median(0,val-1, prof.cycles[i]);
			if ((val % 2) == 1){ // odd
				Q1 = __median(0,val/2 -1, prof.cycles[i]);
				// skip the middle element
				Q3 = __median(val/2 +1,val-1, prof.cycles[i]);
			}
			else{
				Q1 = __median(1,val/2 -2 , prof.cycles[i]);
				// ski the two middle elements
				Q3 = __median(val/2 +1,val-1, prof.cycles[i]);
			}
			len += snprintf(&procfs_buffer[len], 5*80*PROFWIDTH - len, "Profiler Class %d\nMedian (Q2) is %llu; Q1=%llu; Q3=%llu \n", i, Q2, Q1, Q3);

			/* calculate the average */
			average = 0;
			for (j=0; j<val; j++){
				average += ((double) prof.cycles[i][j]) / val;
				//dprintf("%d,%d:%llu\n",i,j,prof.cycles[i][j]);
			}
			len += snprintf(&procfs_buffer[len], 
					5 * 80 * PROFWIDTH - len, 
					"Average is %lf\n\n", average);
		}
	}

	*buffer_location = procfs_buffer;

	/* Return the length */
	return len;
}

static struct proc_dir_entry* profiler_proc_entry=NULL;

/** register to procfs. Automatically calls init_profiler (just in case you forget) */
void __internal_profiler_procfs_open(void){
	if (profiler_proc_entry)
		return;

	profiler_proc_entry = create_proc_read_entry(PROFILER_PROC_NAME, S_IFREG | S_IRUGO, NULL, (read_proc_t*) &__internal_profiler_procfs, NULL);
	profiler_proc_entry->owner = THIS_MODULE;
}

/** unregister from procfs */
void __internal_profiler_procfs_close(void){
	remove_proc_entry(PROFILER_PROC_NAME,NULL);
}
#endif /* __KERNEL__ */

#endif /* PROFILE */
