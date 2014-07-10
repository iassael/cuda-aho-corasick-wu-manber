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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "ixp2000-ue-disas.h"

/* private declarations *****************************************************/
#define BANK_A		0
#define BANK_B		1


/* restricted/unrestricted field decoding ***********************************/
static char *decode_restricted(int indirect_msb, int num, int bank, int is_dest)
{
	static char *index_reg[32] = {
		"*$index", "*n$index", "*$index++", "*n$index++",
		"*$index--", "undef(45)", "undef(46)", "undef(47)",
		"*$$index", "undef(49)", "*$$index++", "undef(4b)",
		"*$$index--", "undef(4d)", "undef(4e)", "undef(4f)",

		"*l$index0", "*l$index0[1]", "*l$index0[2]", "*l$index0[3]",
		"*l$index0[4]", "*l$index0[5]", "*l$index0[6]", "*l$index0[7]",
		"*l$index1", "*l$index1[1]", "*l$index1[2]", "*l$index1[3]",
		"*l$index1[4]", "*l$index1[5]", "*l$index1[6]", "*l$index1[7]",
	};
	char temp[32];

	if (is_dest && num == 0x20) {
		sprintf(temp, "--");
	} else if (num & 0x20) {
		sprintf(temp, "0x%.2x", (indirect_msb ? 0x80 : 0) |
					((num & 0xc0) >> 1) | (num & 0x1f));
	} else if ((num & 0xc0) == 0x00) {
		sprintf(temp, "%c%d", bank==BANK_A?'a':'b', num&0x1f);
	} else if ((num & 0xc0) == 0x40 && !(is_dest && num == 0x41)) {
		sprintf(temp, index_reg[num&0x1f]);
	} else if ((num & 0xc0) == 0x80) {
		sprintf(temp, "$%d", num&0x1f);
	} else if ((num & 0xc0) == 0xc0) {
		sprintf(temp, "$$%d", num&0x1f);
	} else {
		sprintf(temp, "undef(%.2x)", num);
	}

	return strdup(temp);
}

static char *decode_unrestricted(int num, int bank, int is_dest)
{
	char temp[32];

	if ((num & 0x3e0) == 0x000) {
		sprintf(temp, "%c%d", bank==BANK_A?'a':'b', num&0x1f);
	} else if ((num & 0x380) == 0x080) {
		sprintf(temp, "@%c%d", bank==BANK_A?'a':'b', num&0x7f);
	} else if (num == 0x100) {
		sprintf(temp, "*$$index");
	} else if (num == 0x102) {
		sprintf(temp, "*$$index++");
	} else if (num == 0x104) {
		sprintf(temp, "*$$index--");
	} else if (num == 0x140) {
		sprintf(temp, "*$index");
	} else if (num == 0x142) {
		sprintf(temp, "*$index++");
	} else if (num == 0x144) {
		sprintf(temp, "*$index--");
	} else if ((num & 0x3e0) == 0x180) {
		sprintf(temp, "$%d", num & 0x1f);
	} else if ((num & 0x3e0) == 0x1c0) {
		sprintf(temp, "$$%d", num & 0x1f);
	} else if ((num & 0x3f0) == 0x200) {
		sprintf(temp, "*l$index0[%d]", num & 0xf);
	} else if (num == 0x210) {
		sprintf(temp, "*l$index0++");
	} else if (num == 0x211) {
		sprintf(temp, "*l$index0--");
	} else if ((num & 0x3f0) == 0x220) {
		sprintf(temp, "*l$index1[%d]", num & 0xf);
	} else if (num == 0x230) {
		sprintf(temp, "*l$index1++");
	} else if (num == 0x231) {
		sprintf(temp, "*l$index1--");
	} else if (!is_dest && num == 0x241) {
		sprintf(temp, "*n$index");
	} else if (num == 0x243) {
		sprintf(temp, "*n$index++");
	} else if ((num & 0x3e0) == 0x280) {
		sprintf(temp, "n$%d", num & 0x1f);
	} else if (is_dest && num == 0x300) {
		sprintf(temp, "--");
	} else if ((num & 0x300) == 0x300) {
		sprintf(temp, "0x%.2x", num & 0xff);
	} else {
		sprintf(temp, "undef(%.3x)", num);
	}

	return strdup(temp);
}


/* MEM opcode decoding ******************************************************/
#define XFER_TYPE_SRAM	0
#define XFER_TYPE_DRAM	1

static char *decode_nn(int nn, int sig)
{
	char temp[32];

	if (sig == 0 && nn == 3) {
		strcpy(temp, "");
	} else if (nn == 0) {
		sprintf(temp, ", ctx_swap[s%d]", sig);
	} else if (nn == 1) {
		sprintf(temp, ", ctx_swap[s%d], defer[1]", sig);
	} else if (nn == 2) {
		sprintf(temp, ", ctx_swap[s%d], defer[2]", sig);
	} else if (nn == 3) {
		sprintf(temp, ", sig_done[s%d]", sig);
	} else {
		sprintf(temp, ", unknown");
	}

	return strdup(temp);
};

static char *decode_xfer_mem(int num, int contexts_mode, int default_xfer_type)
{
	char temp[32];
	int xfer_type;

	xfer_type = default_xfer_type;
	if (contexts_mode == CONTEXTS_8 && num & 0x10) {
		xfer_type ^= 1;
		num &= 0x0f;
	}

	sprintf(temp, "%s%d", (xfer_type == XFER_TYPE_DRAM) ? "$$" : "$", num);

	return strdup(temp);
}

static char *sram(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	int cmd;
	char *xfer;
	char *a_op;
	char *b_op;
	int ref_cnt;
	char *nn;
	char *indirect_ref;
	char *ignore_data_error;

	cmd = (insn >> 32) & 0xf;
	xfer = decode_xfer_mem((insn >> 20) & 0x1f, contexts_mode, XFER_TYPE_SRAM);
	a_op = decode_restricted(0, insn & 0xff, BANK_A, 0);
	b_op = decode_restricted(0, (insn >> 10) & 0xff, BANK_B, 0);
	ref_cnt = ((insn >> 25) & 7) + 1;
	nn = decode_nn((insn >> 8) & 3, (insn >> 28) & 0xf);
	indirect_ref = ((insn >> 38) & 1) ? ", indirect_ref" : "";
	ignore_data_error = ((insn >> 19) & 1) ? ", ignore_data_error" : "";

	if (cmd == 0 || cmd == 1) {
		sprintf(temp, "sram[%s, %s, %s, %s, %d]%s%s%s",
				(cmd == 0) ? "read" : "write",
				xfer, a_op, b_op, ref_cnt,
				nn, indirect_ref, ignore_data_error);
	} else if (cmd == 2) {
		sprintf(temp, "sram[swap, %s, %s, %s]%s%s",
				xfer, a_op, b_op, nn, indirect_ref);
	} else if (cmd == 3 || cmd == 4 || cmd == 5) {
		sprintf(temp, "sram[%s%s, %s, %s, %s]%s%s",
				((insn >> 18) & 1) ? "test_and_" : "",
				(cmd == 3) ? "set" :
				 (cmd == 4) ? "clr" :
				 (cmd == 5) ? "add" : "unknown",
				xfer, a_op, b_op, nn, indirect_ref);
	} else if (cmd == 7 || cmd == 8) {
		sprintf(temp, "sram[%s%s, %s, %s, %s]%s%s",
				((insn >> 18) & 1) ? "test_and_" : "",
				(cmd == 7) ? "incr" :
				 (cmd == 8) ? "decr" : "unknown",
				((insn >> 18) & 1) ? xfer : "--",
				a_op, b_op, nn, indirect_ref);
	} else if (cmd == 9) {
		sprintf(temp, "sram[%s, %s, %s, %s, %d]%s%s",
				((insn >> 18) & 1) ? "get" : "put",
				xfer, a_op, b_op, ref_cnt,
				nn, indirect_ref);
	} else if (cmd == 10 && ((insn >> 18) & 1) == 0) {
		sprintf(temp, "sram[journal, %s, %s, %s, %d]%s%s",
				xfer, a_op, b_op, ref_cnt,
				nn, indirect_ref);
	} else if (cmd == 10 && ((insn >> 18) & 1) == 1) {
		sprintf(temp, "sram[fast_journal, --, %s, %s]%s%s",
				a_op, b_op, nn, indirect_ref);
	} else if (cmd == 11) {
		sprintf(temp, "sram[dequeue, %s, %s, %s]%s%s",
				xfer, a_op, b_op, nn, indirect_ref);
	} else if (cmd == 12) {
		sprintf(temp, "sram[enqueue%s, %s, %s, %s]%s%s",
				((insn >> 18) & 1) ? "_tail" : "",
				xfer, a_op, b_op, nn, indirect_ref);
	} else if (cmd == 13) {
		sprintf(temp, "sram[csr_%s, %s, %s, %s]%s%s",
				((insn >> 18) & 1) ? "rd" : "wr",
				xfer, a_op, b_op, nn, indirect_ref);
	} else if (cmd == 14) {
		sprintf(temp, "sram[wr_qdesc%s, --, %s, %s]",
				((insn >> 19) & 1) ? "_count" : "",
				a_op, b_op); 
	} else if (cmd == 15 && ((insn >> 18) & 3) == 0) {
		sprintf(temp, "sram[rd_qdesc_other, --, %s, %s]%s%s",
				a_op, b_op, nn, indirect_ref);
	} else if (cmd == 15) {
		int subtype;

		subtype = (insn >> 18) & 3;
		sprintf(temp, "sram[rd_qdesc_%s, %s, %s, %s, %d]%s%s",
				(subtype == 1) ? "tail" :
				 (subtype == 2) ? "head" :
				 "unknown",
				xfer, a_op, b_op, ref_cnt,
				nn, indirect_ref);
	} else {
		sprintf(temp, "sram[unknown]");
	}

	free(nn);
	free(b_op);
	free(a_op);
	free(xfer);

	return strdup(temp);
}

static char *scratch(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	int cmd;
	char *xfer;
	char *a_op;
	char *b_op;
	int ref_cnt;
	char *nn;
	char *indirect_ref;

	cmd = (insn >> 32) & 0xf;
	xfer = decode_xfer_mem((insn >> 20) & 0x1f, contexts_mode, XFER_TYPE_SRAM);
	a_op = decode_restricted(0, insn & 0xff, BANK_A, 0);
	b_op = decode_restricted(0, (insn >> 10) & 0xff, BANK_B, 0);
	ref_cnt = ((insn >> 25) & 7) + 1;
	nn = decode_nn((insn >> 8) & 3, (insn >> 28) & 0xf);
	indirect_ref = ((insn >> 38) & 1) ? ", indirect_ref" : "";

	if (cmd == 0 || cmd == 1) {
		sprintf(temp, "scratch[%s, %s, %s, %s, %d]%s%s",
				(cmd == 0) ? "read" : "write",
				xfer, a_op, b_op, ref_cnt,
				nn, indirect_ref);
	} else if (cmd == 2) {
		sprintf(temp, "scratch[swap, %s, %s, %s]%s%s",
				xfer, a_op, b_op, nn, indirect_ref);
	} else if (cmd == 3 || cmd == 4 || cmd == 5 || cmd == 6) {
		sprintf(temp, "scratch[%s%s, %s, %s, %s]%s%s",
				((insn >> 18) & 1) ? "test_and_" : "",
				(cmd == 3) ? "set" :
				 (cmd == 4) ? "clr" :
				 (cmd == 5) ? "add" :
				 (cmd == 6) ? "sub" : "unknown",
				xfer, a_op, b_op, nn, indirect_ref);
	} else if (cmd == 7 || cmd == 8) {
		sprintf(temp, "scratch[%s%s, %s, %s, %s]%s%s",
				((insn >> 18) & 1) ? "test_and_" : "",
				(cmd == 7) ? "incr" :
				 (cmd == 8) ? "decr" : "unknown",
				((insn >> 18) & 1) ? xfer : "--",
				a_op, b_op, nn, indirect_ref);
	} else if (cmd == 9 || cmd == 10) {
		sprintf(temp, "scratch[%s, %s, %s, %s, %d]%s%s",
				(cmd == 9) ? "get" : "put",
				xfer, a_op, b_op, ref_cnt,
				nn, indirect_ref);
	} else {
		sprintf(temp, "scratch[unknown]");
	}

	free(nn);
	free(b_op);
	free(a_op);
	free(xfer);

	return strdup(temp);
}

static char *dram(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	int cmd;
	char *xfer;
	char *a_op;
	char *b_op;
	int ref_cnt;
	char *nn;
	char *indirect_ref;
	char *ignore_data_error;

	cmd = (insn >> 33) & 1;
	xfer = decode_xfer_mem((insn >> 20) & 0x1f, contexts_mode, XFER_TYPE_DRAM);
	a_op = decode_restricted(0, insn & 0xff, BANK_A, 0);
	b_op = decode_restricted(0, (insn >> 10) & 0xff, BANK_B, 0);
	ref_cnt = ((insn >> 25) & 7) + 1;
	nn = decode_nn((insn >> 8) & 3, (insn >> 28) & 0xf);
	indirect_ref = ((insn >> 38) & 1) ? ", indirect_ref" : "";
	ignore_data_error = ((insn >> 19) & 1) ? ", ignore_data_error" : "";

	if (cmd == 0) {
		sprintf(temp, "dram[%s, %s, %s, %s, %d]%s%s%s",
				((insn >> 32) & 1) ? "write" : "read",
				xfer, a_op, b_op, ref_cnt,
				nn, indirect_ref, ignore_data_error);
	} else if (cmd == 1) {
		sprintf(temp, "dram[%s, --, %s, %s, %d]%s%s%s",
				((insn >> 32) & 1) ? "tbuf_wr" : "rbuf_rd",
				a_op, b_op, ref_cnt,
				nn, indirect_ref, ignore_data_error);
	}

	free(nn);
	free(b_op);
	free(a_op);
	free(xfer);

	return strdup(temp);
}

static char *csr(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	int cmd;
	char *xfer;
	int csr_addr;
	int ref_cnt;
	char *nn;
	char *indirect_ref;

	cmd = (insn >> 33) & 1;
	xfer = decode_xfer_mem((insn >> 20) & 0x1f, contexts_mode, XFER_TYPE_SRAM);
	csr_addr = ((insn >> 2) & 0xff00) | (insn & 0xff);
	ref_cnt = ((insn >> 25) & 7) + 1;
	nn = decode_nn((insn >> 8) & 3, (insn >> 28) & 0xf);
	indirect_ref = ((insn >> 38) & 1) ? ", indirect_ref" : "";

	if (cmd == 0) {
		sprintf(temp, "cap[%s, %s, csr_%.4x, %d]%s%s",
				((insn >> 32) & 1) ? "write" : "read",
				xfer, csr_addr, ref_cnt,
				nn, indirect_ref);
	} else if (cmd == 1) {
		char xfer_data[8];

		if ((insn >> 32) & 1)
			sprintf(xfer_data, "%.4x", (u_int32_t)((insn >> 18) & 0x3fff));
		else
			sprintf(xfer_data, "ALU");

		sprintf(temp, "cap[fast_wr, %s, csr_%.4x]%s",
				xfer_data, csr_addr, indirect_ref);
	}

	free(nn);
	free(xfer);

	return strdup(temp);
}

static char *pci(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	char *xfer;
	char *a_op;
	char *b_op;
	int ref_cnt;
	char *nn;
	char *indirect_ref;

	xfer = decode_xfer_mem((insn >> 20) & 0x1f, contexts_mode, XFER_TYPE_SRAM);
	a_op = decode_restricted(0, insn & 0xff, BANK_A, 0);
	b_op = decode_restricted(0, (insn >> 10) & 0xff, BANK_B, 0);
	ref_cnt = ((insn >> 25) & 7) + 1;
	nn = decode_nn((insn >> 8) & 3, (insn >> 28) & 0xf);
	indirect_ref = ((insn >> 38) & 1) ? ", indirect_ref" : "";

	sprintf(temp, "pci[%s, %s, %s, %s, %d]%s%s",
			((insn >> 32) & 1) ? "write" : "read",
			xfer, a_op, b_op, ref_cnt,
			nn, indirect_ref);

	free(nn);
	free(b_op);
	free(a_op);
	free(xfer);

	return strdup(temp);
}

static char *hash(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	int type;
	char *xfer;
	int ref_cnt;
	char *nn;
	char *indirect_ref;

	type = (insn >> 18) & 3;
	xfer = decode_xfer_mem((insn >> 20) & 0x1f, contexts_mode, XFER_TYPE_SRAM);
	ref_cnt = (insn >> 25) & 3;
	nn = decode_nn((insn >> 8) & 3, (insn >> 28) & 0xf);
	indirect_ref = ((insn >> 38) & 1) ? ", indirect_ref" : "";

	sprintf(temp, "hash_%s[%s, %d]%s%s",
			(type == 0) ? "48" :
			 (type == 1) ? "64" :
			 (type == 2) ? "128" : "unknown",
			xfer, ref_cnt, nn, indirect_ref);

	free(nn);
	free(xfer);

	return strdup(temp);
}

static char *msf(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	int cmd;
	int subtype;
	char *xfer;
	char *a_op;
	char *b_op;
	int ref_cnt;
	char *nn;
	char *indirect_ref;

	cmd = (insn >> 32) & 1;
	subtype = (insn >> 18) & 3;
	xfer = decode_xfer_mem((insn >> 20) & 0x1f, contexts_mode, XFER_TYPE_SRAM);
	a_op = decode_restricted(0, insn & 0xff, BANK_A, 0);
	b_op = decode_restricted(0, (insn >> 10) & 0xff, BANK_B, 0);
	ref_cnt = ((insn >> 25) & 7) + 1;
	nn = decode_nn((insn >> 8) & 3, (insn >> 28) & 0xf);
	indirect_ref = ((insn >> 38) & 1) ? ", indirect_ref" : "";

	if (cmd == 0 || subtype != 1) {
		sprintf(temp, "msf[%s%s, %s, %s, %s, %d]%s%s",
				(cmd == 0) ? "read" : "write",
				(subtype == 0) ? "" :
				 (subtype == 2) ? "64" : "unknown",
				xfer, a_op, b_op, ref_cnt,
				nn, indirect_ref);
	} else if (cmd == 1 && subtype == 1) {
		sprintf(temp, "msf[fast_wr, --, %s, %s]", a_op, b_op);
	}

	free(nn);
	free(b_op);
	free(a_op);
	free(xfer);

	return strdup(temp);
}

static char *cap(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	char *xfer;
	char *a_op;
	char *b_op;
	int ref_cnt;
	char *nn;
	char *indirect_ref;

	xfer = decode_xfer_mem((insn >> 20) & 0x1f, contexts_mode, XFER_TYPE_SRAM);
	a_op = decode_restricted(0, insn & 0xff, BANK_A, 0);
	b_op = decode_restricted(0, (insn >> 10) & 0xff, BANK_B, 0);
	ref_cnt = ((insn >> 25) & 7) + 1;
	nn = decode_nn((insn >> 8) & 3, (insn >> 28) & 0xf);
	indirect_ref = ((insn >> 38) & 1) ? ", indirect_ref" : "";

	sprintf(temp, "cap[%s, %s, %s, %s, %d]%s%s",
			((insn >> 32) & 1) ? "write" : "read",
			xfer, a_op, b_op, ref_cnt,
			nn, indirect_ref);

	free(nn);
	free(b_op);
	free(a_op);
	free(xfer);

	return strdup(temp);
}

/* non-MEM opcode decoding **************************************************/
static char *decode_alu_shf(int type, int amount)
{
	char temp[32];

	if (type == 0) {
		sprintf(temp, "<<rot%d", (32-amount) & 31);
	} else if (type == 1 && amount == 0) {
		sprintf(temp, ">>indirect");
	} else if (type == 1 && amount != 0) {
		sprintf(temp, ">>%d", amount);
	} else if (type == 2 && amount == 0) {
		sprintf(temp, "<<indirect");
	} else if (type == 2 && amount != 0) {
		sprintf(temp, "<<%d", (32-amount) & 31);
	} else if (type == 3) {
		sprintf(temp, ">>%d", amount);
	} else {
		sprintf(temp, "unknown");
	}

	return strdup(temp);
}

static char *decode_defer(int defer)
{
	if (defer == 0) {
		return "";
	} else if (defer == 1) {
		return ", defer[1]";
	} else if (defer == 2) {
		return ", defer[2]";
	} else if (defer == 3) {
		return ", defer[3]";
	}

	return "unknown";
}

static char *alu_shf(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	char *dest;
	char *a_op;
	int alu_op;
	char *b_op;
	char *b_op_shf;

	dest = decode_restricted(0, (insn >> 20) & 0xff, (insn >> 36) & 1, 1);
	a_op = decode_restricted((insn >> 18) & 1, insn & 0xff, BANK_A, 0);
	alu_op = (insn >> 33) & 7;
	b_op = decode_restricted((insn >> 18) & 1, (insn >> 10) & 0xff, BANK_B, 0);
	b_op_shf = decode_alu_shf((insn >> 8) & 3, (insn >> 28) & 0x1f);

	if ((insn >> 19) & 1) {
		char *t;

		t = a_op;
		a_op = b_op;
		b_op = t;
	}

	if (alu_op == 0 && ((insn >> 8) & 3) == 3) {
		sprintf(temp, "dbl_shf[%s, %s, %s, %s]",
				dest, a_op, b_op, b_op_shf);
	} else if (alu_op >= 0 && alu_op <= 5) {
		sprintf(temp, "alu_shf[%s, %s, %s, %s, %s]",
				dest, a_op,
				(alu_op == 0) ? "B" :
				 (alu_op == 1) ? "~B" :
				 (alu_op == 2) ? "AND" :
				 (alu_op == 3) ? "~AND" :
				 (alu_op == 4) ? "AND~" :
				 (alu_op == 5) ? "OR" :
				 (alu_op == 6) ? "ASR" :
				 (alu_op == 7) ? "BYTE_ALIGN" : "unknown",
				b_op, b_op_shf);
	} else if (alu_op == 6) {
		sprintf(temp, "asr[%s, %s, %s]", dest, b_op, b_op_shf);
	} else if (alu_op == 7) {
		int subtype;

		subtype = (insn >> 8) & 3;
		sprintf(temp, "byte_align_%s[%s, %s]",
				(subtype == 1) ? "be" :
				 (subtype == 2) ? "le" : "unknown",
				dest, a_op);
	}

	free(b_op_shf);
	free(b_op);
	free(a_op);
	free(dest);

	return strdup(temp);
}

static char *alu(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	char *dest;
	char *destn;
	char *a_op;
	int alu_op;
	char *b_op;
	char *b_op_dest;

	dest = decode_unrestricted((insn >> 20) & 0x3ff, (insn >> 36) & 1, 1);
	destn = decode_unrestricted((insn >> 20) & 0x3ff, (insn >> 36) & 1, 0);
	a_op = decode_unrestricted(insn & 0x3ff, BANK_A, 0);
	alu_op = (insn >> 31) & 0x1f;
	b_op = decode_unrestricted((insn >> 10) & 0x3ff, BANK_B, 0);
	b_op_dest = decode_unrestricted((insn >> 10) & 0x3ff, BANK_B, 1);

	if ((insn >> 30) & 1) {
		char *t;

		t = a_op;
		a_op = b_op;
		b_op = t;
	}

	if (alu_op == 3) {
		sprintf(temp, "pop_count3[%s, %s]", dest, b_op);
	} else if (alu_op == 6 || alu_op == 7) {
		sprintf(temp, "pop_count%d[%s]", alu_op - 5, b_op);
	} else if (alu_op == 11) {
		sprintf(temp, "cam_clear");
	} else if (alu_op == 14) {
		sprintf(temp, "ffs[%s, %s]", dest, b_op);
	} else if (alu_op == 15) {
		sprintf(temp, "cam_read_tag[%s, %s]", dest, b_op);
	} else if (alu_op == 18) {
		int crc_type;
		int crc_bytes;

		crc_type = (insn >> 15) & 7;
		crc_bytes = (insn >> 10) & 7;
		sprintf(temp, "crc_%s[%s, %s, %s], %s%s",
				((insn >> 13) & 1) ? "le" : "be",
				(crc_type == 0) ? "none" :
				 (crc_type == 2) ? "crc_ccitt" :
				 (crc_type == 4) ? "crc_32" :
				 (crc_type == 5) ? "crc_iscsi" :
				 (crc_type == 6) ? "crc_10" :
				 (crc_type == 7) ? "crc_5" : "unknown",
				dest, b_op,
				(crc_bytes == 0) ? "crc_bytes_0_3" :
				 (crc_bytes == 1) ? "crc_bytes_1_3" :
				 (crc_bytes == 2) ? "crc_bytes_2_3" :
				 (crc_bytes == 3) ? "crc_bytes_3" :
				 (crc_bytes == 4) ? "crc_bytes_0_2" :
				 (crc_bytes == 5) ? "crc_bytes_0_1" :
				 (crc_bytes == 6) ? "crc_bytes_0" :
				 "unknown",
				((insn >> 14) & 1) ? ", bit_swap" : "");
	} else if (alu_op == 19) {
		sprintf(temp, "cam_write[%s, %s, %s]",
				b_op_dest, a_op, destn);
	} else if (alu_op == 23) {
		char opt_tok[16];
		int right_op;

		strcpy(opt_tok, "");
		if ((insn >> 30) & 1)
			right_op = insn & 0x3ff;
		else
			right_op = (insn >> 10) & 0x3ff;

		if (right_op & 1) {
			sprintf(opt_tok, ", lm_addr0[%d]", (right_op >> 2) & 3);
		} else if (right_op & 2) {
			sprintf(opt_tok, ", lm_addr1[%d]", (right_op >> 2) & 3);
		}

		sprintf(temp, "cam_lookup[%s, %s]%s", dest, a_op, opt_tok);
	} else if (alu_op == 27) {
		sprintf(temp, "cam_write_state[%s, %s]", b_op_dest, destn);
	} else if (alu_op == 31) {
		sprintf(temp, "cam_read_state[%s, %s]", dest, b_op);
	} else {
		sprintf(temp, "alu[%s, %s, %s, %s]",
				dest, a_op,
				(alu_op == 0) ? "B" :
				 (alu_op == 1) ? "+" :
				 (alu_op == 4) ? "~B" :
				 (alu_op == 5) ? "+16" :
				 (alu_op == 8) ? "AND" :
				 (alu_op == 9) ? "+8" :
				 (alu_op == 12) ? "~AND" :
				 (alu_op == 13) ? "-CARRY" :
				 (alu_op == 16) ? "AND~" :
				 (alu_op == 17) ? "+CARRY" :
				 (alu_op == 20) ? "OR" :
				 (alu_op == 21) ? "-" :
				 (alu_op == 24) ? "XOR" :
				 (alu_op == 25) ? "B-A" : "unknown",
				b_op);
	}

	free(b_op_dest);
	free(b_op);
	free(a_op);
	free(destn);
	free(dest);

	return strdup(temp);
}

static char *ld_field(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	char *dest_reg;
	char byte_enables[8];
	int indirect_msb;
	char *src_op;
	char *op_shf_cntl;
	char *opt_tok;

	indirect_msb = (insn >> 18) & 1;
	if ((insn >> 19) & 1) {
		dest_reg = decode_restricted(indirect_msb,
						(insn >> 10) & 0xff,
						BANK_B, 1);
		src_op = decode_restricted(indirect_msb, insn & 0xff,
						BANK_A, 0);
	} else {
		dest_reg = decode_restricted(indirect_msb, insn & 0xff,
						BANK_A, 1);
		src_op = decode_restricted(indirect_msb, (insn >> 10) & 0xff,
						BANK_B, 0);
	}
	sprintf(byte_enables, "%c%c%c%c",
		((insn >> 27) & 1) ? '1' : '0',
		((insn >> 26) & 1) ? '1' : '0',
		((insn >> 25) & 1) ? '1' : '0',
		((insn >> 24) & 1) ? '1' : '0');
	op_shf_cntl = decode_alu_shf((insn >> 8) & 3, (insn >> 28) & 0x1f);
	opt_tok = ((insn >> 34) & 1) ? ", load_cc" : "";

	sprintf(temp, "ld_field%s[%s, %s, %s, %s]%s",
			((insn >> 20) & 1) ? "_w_clr" : "",
			dest_reg, byte_enables, src_op, op_shf_cntl,
			opt_tok);

	free(op_shf_cntl);
	free(src_op);
	free(dest_reg);

	return strdup(temp);
}

static char *branch_byte(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	int indirect_msb;
	char *reg;
	int byte_no;
	char *byte_compare_value;
	char label[16];
	char *opt_tok;

	indirect_msb = (insn >> 18) & 1;
	if ((insn >> 15) & 1) {
		reg = decode_restricted(indirect_msb, insn & 0x3ff, BANK_A, 1);
		byte_compare_value = decode_restricted(indirect_msb,
							(insn >> 10) & 0x3ff,
							BANK_B, 0);
	} else {
		reg = decode_restricted(indirect_msb, (insn >> 10) & 0x3ff,
					BANK_B, 1);
		byte_compare_value = decode_restricted(indirect_msb,
							insn & 0x3ff,
							BANK_A, 0);
	}
	byte_no = (insn >> 8) & 3;
	sprintf(label, "l%d#", (u_int32_t)((insn >> 22) & 0x1fff));
	opt_tok = decode_defer((insn >> 20) & 3);

	sprintf(temp, "br%s=byte[%s, %d, %s, %s]%s",
			((insn >> 19) & 1) ? "" : "!",
			reg, byte_no, byte_compare_value, label,
			opt_tok);

	free(byte_compare_value);
	free(reg);

	return strdup(temp);
}

static char *branch_bit(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	char *reg;
	int bit_position;
	char label[16];
	char *opt_tok;

	if ((insn >> 15) & 1) {
		reg = decode_restricted(0, insn & 0x3ff, BANK_A, 1);
		bit_position = ((insn >> 10) + 31) & 0x1f;
	} else {
		reg = decode_restricted(0, (insn >> 10) & 0x3ff, BANK_B, 1);
		bit_position = (insn + 31) & 0x1f;
	}
	sprintf(label, "l%d#", (u_int32_t)((insn >> 22) & 0x1fff));
	opt_tok = decode_defer((insn >> 20) & 3);

	sprintf(temp, "br_%s[%s, %d, %s]%s",
			((insn >> 18) & 1) ? "bset" : "bclr",
			reg, bit_position, label, opt_tok);

	free(reg);

	return strdup(temp);
}

static char *branch(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	char label[16];
	int parameter;
	int pipe_stage;
	int condition;
	char *opt_tok;

	sprintf(label, "l%d#", (u_int32_t)((insn >> 22) & 0x1fff));
	parameter = (insn >> 14) & 0xf;
	pipe_stage = ((insn >> 8) & 3) + 1;
	condition = insn & 0x1f;
	opt_tok = decode_defer((insn >> 20) & 3);

	if (condition >= 0 && condition <= 11) {
		sprintf(temp, "b%s[%s], stage[%d]%s",
			(condition == 0) ? "eq" :
			 (condition == 1) ? "ne" :
			 (condition == 2) ? "mi" :
			 (condition == 3) ? "pl" :
			 (condition == 4) ? "cs/bhs" :
			 (condition == 5) ? "cc/blo" :
			 (condition == 6) ? "vs" :
			 (condition == 7) ? "vc" :
			 (condition == 8) ? "gt" :
			 (condition == 9) ? "lt" :
			 (condition == 10) ? "le" :
			 (condition == 11) ? "gt" : "unknown",
			label, pipe_stage, opt_tok);
	} else if (condition >= 16 && condition <= 17) {
		sprintf(temp, "br%s=ctx[%d, %s]%s",
			(condition == 16) ? "" : "!",
			parameter, label, opt_tok);
	} else if (condition >= 18 && condition <= 19) {
		sprintf(temp, "br_%ssignal[s%d, %s]%s",
			(condition == 18) ? "" : "!",
			parameter, label, opt_tok);
	} else if (condition >= 20 && condition <= 21) {
		sprintf(temp, "br_%sinp_state[%s, %s]%s",
			(condition == 20) ? "" : "!",
			(parameter == 0) ? "NN_EMPTY" :
			 (parameter == 1) ? "NN_FULL" :
			 (parameter == 2) ? "SCR_RING0_STATUS" :
			 (parameter == 3) ? "SCR_RING1_STATUS" :
			 (parameter == 4) ? "SCR_RING2_STATUS" :
			 (parameter == 5) ? "SCR_RING3_STATUS" :
			 (parameter == 6) ? "SCR_RING4_STATUS" :
			 (parameter == 7) ? "SCR_RING5_STATUS" :
			 (parameter == 8) ? "SCR_RING6_STATUS" :
			 (parameter == 9) ? "SCR_RING7_STATUS" :
			 (parameter == 10) ? "SCR_RING8_STATUS" :
			 (parameter == 11) ? "SCR_RING9_STATUS" :
			 (parameter == 12) ? "SCR_RING10_STATUS" :
			 (parameter == 13) ? "SCR_RING11_STATUS" :
			 (parameter == 14) ? "FCI_NOT_EMPTY" :
			 (parameter == 15) ? "FCI_FULL" : "unknown",
			label, opt_tok);
	} else if (condition == 24) {
		sprintf(temp, "br[%s]%s", label, opt_tok);
	} else {
		sprintf(temp, "br[unknown]");
	}

	return strdup(temp);
}

static char *context(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	char signal_list[128];
	char opt_tok[128];
	int defer;
	int i;

	memset(signal_list, 0, 4);
	for (i=1;i<=15;i++) {
		if ((insn >> i) & 1) {
			char l[8];

			sprintf(l, ", s%d", i);
			strcat(signal_list, l);
		}
	}
	if (insn & 1)
		strcat(signal_list, ", voluntary");
	if ((insn & 0x1ffff) == 0x10000)
		strcat(signal_list, ", kill");
	if ((insn >> 17) & 1)
		strcat(signal_list, ", bpt");
	if ((insn >> 19) & 1)
		strcat(signal_list, ", --");

	strcpy(opt_tok, "");
	if ((insn >> 16) & 1) {
		strcat(opt_tok, ", any");
	}
	defer = (insn >> 20) & 3;
	if (defer == 1) {
		strcat(opt_tok, ", defer[1]");
	} else if (defer == 2) {
		strcat(opt_tok, ", defer[2]");
	}
	if ((insn >> 18) & 1) {
		char label[16];

		sprintf(label, ", br[l%d#]", (u_int32_t)((insn >> 22) & 0x1fff));
		strcat(opt_tok, label);
	}

	sprintf(temp, "ctx_arb[%s]%s", signal_list + 2, opt_tok);

	return strdup(temp);
}

static char *branch_alu(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	char *src;
	int address;

	if (((insn >> 18) & 3) == 3) {
		address = (insn >> 10) & 0xff;
		src = decode_unrestricted(insn & 0x3ff, BANK_A, 0);
	} else {
		address = insn & 0xff;
		src = decode_unrestricted((insn >> 10) & 0x3ff, BANK_B, 0);
	}
	address |= ((insn >> 22) & 0x1f) << 8;

	sprintf(temp, "jump[%s, l%d#]%s",
			src, address, decode_defer((insn >> 20) & 3));

	free(src);

	return strdup(temp);
}

static char *immed(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	int type;
	char *dest;
	u_int32_t immed;

	if (((insn >> 18) & 3) == 3) {
		immed = (insn >> 10) & 0xff;
		dest = decode_unrestricted(insn & 0x3ff, BANK_A, 1);
	} else {
		immed = insn & 0xff;
		dest = decode_unrestricted((insn >> 10) & 0x3ff, BANK_B, 1);
	}

	type = (insn >> 29) & 3;
	if (type == 0) {
		char shift[16];

		immed |= ((insn >> 20) & 0xff) << 8;
		if ((insn >> 31) & 1)
			immed = ~immed;
		sprintf(shift, "<<%d", 8 * (u_int32_t)((insn >> 33) & 3));

		sprintf(temp, "immed[%s, 0x%.8x, %s]", dest, immed, shift);
	} else if (type == 1) {
		sprintf(temp, "immed_b%d[%s, 0x%.2x]",
				(u_int32_t)((insn >> 33) & 3), dest, immed);
	} else if (type == 2) {
		immed |= ((insn >> 20) & 0xff) << 8;
		sprintf(temp, "immed_w%d[%s, 0x%.4x]",
				(u_int32_t)((insn >> 34) & 1), dest, immed);
	} else {
		sprintf(temp, "immed[unknown]");
	}

	free(dest);

	return strdup(temp);
}

static char *multiply(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	char *type;

	switch ((insn >> 31) & 3) {
	case 0:	type = "start"; break;
	case 1:	type = "24x8"; break;
	case 2:	type = "16x16"; break;
	case 3:	type = "32x32"; break;
	}

	if ((insn >> 22) & 1) {
		char *dest;

		dest = decode_unrestricted(insn & 0x3ff, (insn >> 23) & 1, 1);
		sprintf(temp, "mul_step[%s, --], %s_last%s",
				dest, type, (insn >> 20) ? "2" : "");
		free(dest);
	} else {
		char *a_op;
		char *b_op;
		int subtype;

		a_op = decode_unrestricted(insn & 0x3ff, BANK_A, 0);
		b_op = decode_unrestricted((insn >> 10) & 0x3ff, BANK_B, 0);
		if ((insn >> 30) & 1) {
			char *t;

			t = a_op;
			a_op = b_op;
			b_op = t;
		}
		subtype = (insn >> 20) & 7;

		sprintf(temp, "mul_step[%s, %s], %s%s",
				a_op, b_op, type,
				((insn >> 31) & 3) ?
				 (subtype == 0) ? "_step1" :
				 (subtype == 1) ? "_step2" :
				 (subtype == 2) ? "_step3" :
				 (subtype == 3) ? "_step4" : "unknown" : "");

		free(b_op);
		free(a_op);
	}

	return strdup(temp);
}

static char *local_csr(u_int64_t insn, int contexts_mode)
{
	char temp[128];
	int csr_addr;

	csr_addr = (insn >> 22) & 0x7ff;
	if (((insn >> 21) & 1) == 0) {
		sprintf(temp, "local_csr_rd[local_csr_%.3x]", csr_addr);
	} else {
		char *src;

		if (((insn >> 8) & 3) == 3)
			src = decode_unrestricted((insn >> 10) & 0x3ff, BANK_B, 0);
		else
			src = decode_unrestricted(insn & 0x3ff, BANK_A, 0);

		sprintf(temp, "local_csr_wr[local_csr_%.3x, %s]", csr_addr, src);

		free(src);
	}

	return strdup(temp);
}


/* frontend *****************************************************************/
char *ixp2000_ue_disassemble(u_int64_t insn, int contexts_mode)
{
	u_int8_t msb;

	msb = (insn >> 32) & 0xff;
	if ((msb & 0xb0) == 0x00)		/* 0-00 ---- */
		return sram(insn, contexts_mode);
	else if ((msb & 0xb0) == 0x10)		/* 0-01 ---- */
		return scratch(insn, contexts_mode);
	else if ((msb & 0xbc) == 0x28)		/* 0-10 10-- */
		return dram(insn, contexts_mode);		
	else if ((msb & 0xbc) == 0x30)		/* 0-11 00-- */
		return csr(insn, contexts_mode);
	else if ((msb & 0xbe) == 0x34)		/* 0-11 010- */
		return pci(insn, contexts_mode);
	else if ((msb & 0xbe) == 0x36)		/* 0-11 011- */
		return hash(insn, contexts_mode);
	else if ((msb & 0xbe) == 0x3c)		/* 0-11 110- */
		return msf(insn, contexts_mode);
	else if ((msb & 0xbe) == 0x3e)		/* 0-11 111- */
		return cap(insn, contexts_mode);
	else if ((msb & 0xe0) == 0x80)		/* 100- ---- */
		return alu_shf(insn, contexts_mode);
	else if ((msb & 0xe0) == 0xa0)		/* 101- ---- */
		return alu(insn, contexts_mode);
	else if ((msb & 0xf8) == 0xc0)		/* 1100 0--- */
		return ld_field(insn, contexts_mode);
	else if ((msb & 0xf8) == 0xc8)		/* 1100 1--- */
		return branch_byte(insn, contexts_mode);
	else if ((msb & 0xf8) == 0xd0)		/* 1101 0--- */
		return branch_bit(insn, contexts_mode);
	else if ((msb & 0xf8) == 0xd8)		/* 1101 1--- */
		return branch(insn, contexts_mode);
	else if ((msb & 0xf8) == 0xe0)		/* 1110 0--- */
		return context(insn, contexts_mode);
	else if ((msb & 0xf8) == 0xe8)		/* 1110 1--- */
		return branch_alu(insn, contexts_mode);
	else if ((msb & 0xf8) == 0xf0)		/* 1111 0--- */
		return immed(insn, contexts_mode);
	else if ((msb & 0xfe) == 0xf8)		/* 1111 100- */
		return multiply(insn, contexts_mode);
	else if ((msb & 0xfe) == 0xfc)		/* 1111 110- */
		return local_csr(insn, contexts_mode);

	return NULL;
}


#ifdef TEST
int main(int argc, char *argv[])
{
	int i;

	for (i=1;i<argc;i++) {
		u_int64_t insn;

		if (sscanf(argv[i], "%Lx", &insn) == 1) {
			char *d;

			d = ixp2000_ue_disassemble(insn, CONTEXTS_8);
			printf("%.10Lx = %s\n", insn, d);
			free(d);
		}
	}

	return 0;
}
#endif
