/*
 * Generic library functions for the MSF (Media and Switch Fabric
 * unit) and microengines found on the Intel IXP2000 series of network
 * processors.
 *
 * Stub functions to make it work from userspace.
 *
 * Copyright (C) 2004, 2005 Lennert Buytenhek <buytenh@wantstofly.org>
 * Dedicated to Marija Kulikova.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include "compat.h"

#define dprintf(...)

static int dev_mem_fd;
void *IXP2000_GLOBAL_REG_VIRT_BASE;
void *IXP2000_MSF_VIRT_BASE;
void *IXP2000_RBUF_TBUF_VIRT_BASE;
void *IXP2000_UENGINE_CSR_VIRT_BASE;
void *IXP2000_INT_CONTROLLER_VIRT_BASE;
u32 ixp2000_uengine_mask;

static void ixp2000_map(void) __attribute__((constructor));
static void ixp2000_map(void)
{
	u32 product_id;

	dev_mem_fd = open("/dev/mem", O_RDWR | O_SYNC);
	if (dev_mem_fd < 0) {
		perror("open(\"/dev/mem\")");
		exit(-1);
	}

	IXP2000_GLOBAL_REG_VIRT_BASE = ioremap_nocache(0xc0004000, 4096);
	IXP2000_MSF_VIRT_BASE = ioremap_nocache(0xc8000000, 8192);
	IXP2000_RBUF_TBUF_VIRT_BASE = ioremap_nocache(0xc8002000, 8192);
	IXP2000_UENGINE_CSR_VIRT_BASE = ioremap_nocache(0xc0018000, 32768);
	IXP2000_INT_CONTROLLER_VIRT_BASE = ioremap_nocache(0xd6000000, 4096);

	// @@@ we should check that we're really on an ixp2000
	product_id = *IXP2000_PRODUCT_ID;

	switch ((product_id >> 8) & 0x1fff) {
	case 0:
		dprintf("detected IXP2800 rev %c%x\n",
			'A' + ((product_id >> 4) & 0xf), product_id & 0xf);
		ixp2000_uengine_mask = 0x00ff00ff;
		break;

	case 1:
		dprintf("detected IXP2850 rev %c%x\n",
			'A' + ((product_id >> 4) & 0xf), product_id & 0xf);
		ixp2000_uengine_mask = 0x00ff00ff;
		break;

	case 2:
		dprintf("detected IXP2400 rev %c%x\n",
			'A' + ((product_id >> 4) & 0xf), product_id & 0xf);
		ixp2000_uengine_mask = 0x000f000f;
		break;

	default:
		fprintf(stderr, "unknown ixp2000 model (%.8x)\n", product_id);
		ixp2000_uengine_mask = 0;
		break;
	}
}

static void ixp2000_unmap(void) __attribute__((destructor));
static void ixp2000_unmap(void)
{
	if (dev_mem_fd >= 0) {
		iounmap_length(IXP2000_GLOBAL_REG_VIRT_BASE, 4096);
		iounmap_length(IXP2000_MSF_VIRT_BASE, 8192);
		iounmap_length(IXP2000_RBUF_TBUF_VIRT_BASE, 8192);
		iounmap_length(IXP2000_UENGINE_CSR_VIRT_BASE, 32768);
		iounmap_length(IXP2000_INT_CONTROLLER_VIRT_BASE, 4096);
		close(dev_mem_fd);
	}
}

void *ioremap_nocache(unsigned long phys, unsigned long size)
{
	void *x;

	x = mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, dev_mem_fd, phys);
	if (x == MAP_FAILED) {
		perror("mmap");
		exit(-1);
	}

	return x;
}

void iounmap_length(volatile void *virt, unsigned long size)
{
	munmap((void *)virt, size);
}

void udelay(unsigned long usecs)
{
	usleep(usecs);
}
