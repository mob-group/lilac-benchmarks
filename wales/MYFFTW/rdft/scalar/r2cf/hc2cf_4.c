/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Fri Jan 27 16:12:26 EST 2017 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2c.native -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 4 -dit -name hc2cf_4 -include hc2cf.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 16 additions, 6 multiplications, 6 fused multiply/add),
 * 31 stack variables, 0 constants, and 16 memory accesses
 */
#include "hc2cf.h"

static void hc2cf_4(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 6); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 6, MAKE_VOLATILE_STRIDE(16, rs)) {
	       E To, Te, Tm, T8, Tw, Ty, Tq, Tk;
	       {
		    E T1, Tv, Tu, T7, Tg, Tj, Tf, Ti, Tp, Th;
		    T1 = Rp[0];
		    Tv = Rm[0];
		    {
			 E T3, T6, T2, T5;
			 T3 = Rp[WS(rs, 1)];
			 T6 = Rm[WS(rs, 1)];
			 T2 = W[2];
			 T5 = W[3];
			 {
			      E Ta, Td, Tc, Tn, Tb, Tt, T4, T9;
			      Ta = Ip[0];
			      Td = Im[0];
			      Tt = T2 * T6;
			      T4 = T2 * T3;
			      T9 = W[0];
			      Tc = W[1];
			      Tu = FNMS(T5, T3, Tt);
			      T7 = FMA(T5, T6, T4);
			      Tn = T9 * Td;
			      Tb = T9 * Ta;
			      Tg = Ip[WS(rs, 1)];
			      Tj = Im[WS(rs, 1)];
			      To = FNMS(Tc, Ta, Tn);
			      Te = FMA(Tc, Td, Tb);
			      Tf = W[4];
			      Ti = W[5];
			 }
		    }
		    Tm = T1 - T7;
		    T8 = T1 + T7;
		    Tw = Tu + Tv;
		    Ty = Tv - Tu;
		    Tp = Tf * Tj;
		    Th = Tf * Tg;
		    Tq = FNMS(Ti, Tg, Tp);
		    Tk = FMA(Ti, Tj, Th);
	       }
	       {
		    E Ts, Tr, Tl, Tx;
		    Ts = To + Tq;
		    Tr = To - Tq;
		    Tl = Te + Tk;
		    Tx = Tk - Te;
		    Rp[WS(rs, 1)] = Tm + Tr;
		    Rm[0] = Tm - Tr;
		    Ip[0] = Ts + Tw;
		    Im[WS(rs, 1)] = Ts - Tw;
		    Ip[WS(rs, 1)] = Tx + Ty;
		    Im[0] = Tx - Ty;
		    Rp[0] = T8 + Tl;
		    Rm[WS(rs, 1)] = T8 - Tl;
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 4},
     {TW_NEXT, 1, 0}
};

static const hc2c_desc desc = { 4, "hc2cf_4", twinstr, &GENUS, {16, 6, 6, 0} };

void X(codelet_hc2cf_4) (planner *p) {
     X(khc2c_register) (p, hc2cf_4, &desc, HC2C_VIA_RDFT);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2c.native -compact -variables 4 -pipeline-latency 4 -n 4 -dit -name hc2cf_4 -include hc2cf.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 16 additions, 6 multiplications, 6 fused multiply/add),
 * 13 stack variables, 0 constants, and 16 memory accesses
 */
#include "hc2cf.h"

static void hc2cf_4(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 6); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 6, MAKE_VOLATILE_STRIDE(16, rs)) {
	       E T1, Tp, T6, To, Tc, Tk, Th, Tl;
	       T1 = Rp[0];
	       Tp = Rm[0];
	       {
		    E T3, T5, T2, T4;
		    T3 = Rp[WS(rs, 1)];
		    T5 = Rm[WS(rs, 1)];
		    T2 = W[2];
		    T4 = W[3];
		    T6 = FMA(T2, T3, T4 * T5);
		    To = FNMS(T4, T3, T2 * T5);
	       }
	       {
		    E T9, Tb, T8, Ta;
		    T9 = Ip[0];
		    Tb = Im[0];
		    T8 = W[0];
		    Ta = W[1];
		    Tc = FMA(T8, T9, Ta * Tb);
		    Tk = FNMS(Ta, T9, T8 * Tb);
	       }
	       {
		    E Te, Tg, Td, Tf;
		    Te = Ip[WS(rs, 1)];
		    Tg = Im[WS(rs, 1)];
		    Td = W[4];
		    Tf = W[5];
		    Th = FMA(Td, Te, Tf * Tg);
		    Tl = FNMS(Tf, Te, Td * Tg);
	       }
	       {
		    E T7, Ti, Tn, Tq;
		    T7 = T1 + T6;
		    Ti = Tc + Th;
		    Rm[WS(rs, 1)] = T7 - Ti;
		    Rp[0] = T7 + Ti;
		    Tn = Tk + Tl;
		    Tq = To + Tp;
		    Im[WS(rs, 1)] = Tn - Tq;
		    Ip[0] = Tn + Tq;
	       }
	       {
		    E Tj, Tm, Tr, Ts;
		    Tj = T1 - T6;
		    Tm = Tk - Tl;
		    Rm[0] = Tj - Tm;
		    Rp[WS(rs, 1)] = Tj + Tm;
		    Tr = Th - Tc;
		    Ts = Tp - To;
		    Im[0] = Tr - Ts;
		    Ip[WS(rs, 1)] = Tr + Ts;
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 4},
     {TW_NEXT, 1, 0}
};

static const hc2c_desc desc = { 4, "hc2cf_4", twinstr, &GENUS, {16, 6, 6, 0} };

void X(codelet_hc2cf_4) (planner *p) {
     X(khc2c_register) (p, hc2cf_4, &desc, HC2C_VIA_RDFT);
}
#endif				/* HAVE_FMA */
