/*
 * Disassembler for the IXP2000 microengine (MEv2) instruction format.
 *
 * Copyright (C) 2005 Lennert Buytenhek <buytenh@wantstofly.org>
 * Dedicated to Marija Kulikova.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 */

#ifndef __IXP2000_UE_DISAS_H
#define __IXP2000_UE_DISAS_H

#define CONTEXTS_4	4
#define CONTEXTS_8	8

char *ixp2000_ue_disassemble(u_int64_t insn, int contexts_mode);


#endif
