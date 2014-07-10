#ifndef __COMPAT_H
#define __COMPAT_H

#ifndef __KERNEL__
#include <string.h>
#include <sys/types.h>

typedef u_int8_t u8;
typedef u_int32_t u32;
typedef u_int64_t u64;

extern void *IXP2000_GLOBAL_REG_VIRT_BASE;
extern void *IXP2000_MSF_VIRT_BASE;
extern void *IXP2000_RBUF_TBUF_VIRT_BASE;
extern void *IXP2000_UENGINE_CSR_VIRT_BASE;
extern void *IXP2000_INT_CONTROLLER_VIRT_BASE;
extern u32 ixp2000_uengine_mask;

#define IXP2000_PRODUCT_ID	((volatile u32 *)(IXP2000_GLOBAL_REG_VIRT_BASE + 0x0a00))
#define IXP2000_MISC_CONTROL	((volatile u32 *)(IXP2000_GLOBAL_REG_VIRT_BASE + 0x0a04))
#define IXP2000_MSF_CLK_CNTRL	((volatile u32 *)(IXP2000_GLOBAL_REG_VIRT_BASE + 0x0a08))
#define IXP2000_RESET0		((volatile u32 *)(IXP2000_GLOBAL_REG_VIRT_BASE + 0x0a0c))
#define IXP2000_RESET1		((volatile u32 *)(IXP2000_GLOBAL_REG_VIRT_BASE + 0x0a10))
#define IXP2000_CLOCK_CONTROL	((volatile u32 *)(IXP2000_GLOBAL_REG_VIRT_BASE + 0x0a14))
#define IXP2000_STRAP_OPTIONS	((volatile u32 *)(IXP2000_GLOBAL_REG_VIRT_BASE + 0x0a18))

void *ioremap_nocache(unsigned long phys, unsigned long size);
void iounmap_length(volatile void *virt, unsigned long size);
void udelay(unsigned long usecs);

static inline u32 hweight32(u32 w)
{
	u32 res;

	res = (w & 0x55555555) + ((w >> 1) & 0x55555555);
	res = (res & 0x33333333) + ((res >> 2) & 0x33333333);
	res = (res & 0x0F0F0F0F) + ((res >> 4) & 0x0F0F0F0F);
	res = (res & 0x00FF00FF) + ((res >> 8) & 0x00FF00FF);
	res = (res & 0x0000FFFF) + ((res >> 16) & 0x0000FFFF);

	return res;
}

static inline unsigned int ixdp2x00_master_npu(void)
{
        return !!(*IXP2000_STRAP_OPTIONS & 4);
}
#else
#include <asm/arch-ixp2000/ixp2000-regs.h>
#include <asm/io.h>
#endif


#endif
