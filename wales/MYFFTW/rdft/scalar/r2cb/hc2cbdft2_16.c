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
/* Generated on Fri Jan 27 16:13:31 EST 2017 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2cdft.native -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 16 -dif -name hc2cbdft2_16 -include hc2cb.h */

/*
 * This function contains 206 FP additions, 100 FP multiplications,
 * (or, 136 additions, 30 multiplications, 70 fused multiply/add),
 * 97 stack variables, 3 constants, and 64 memory accesses
 */
#include "hc2cb.h"

static void hc2cbdft2_16(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 30); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 30, MAKE_VOLATILE_STRIDE(64, rs)) {
	       E T3w, T3z, T2Y, T3D, T3x, T3m, T3u, T3C, T3y, T3o, T3k, T3E, T3A;
	       {
		    E T20, Tf, T3Q, T32, T3V, T3f, T2a, TN, T2f, T1m, T3G, T2G, T3L, T2T, T26;
		    E T1F, T3M, T2N, T3H, T2W, T25, Tu, T1n, T1o, T3R, T3i, T2g, T1a, T21, T1y;
		    E T3W, T39;
		    {
			 E T2R, T1B, T2S, T1E;
			 {
			      E T1e, T3, T1C, TA, Tx, T6, T1D, T1h, Td, T1A, TL, T1k, Ta, TC, TF;
			      E T1z;
			      {
				   E T4, T5, T1f, T1g;
				   {
					E T1, T2, Ty, Tz;
					T1 = Rp[0];
					T2 = Rm[WS(rs, 7)];
					Ty = Ip[0];
					Tz = Im[WS(rs, 7)];
					T4 = Rp[WS(rs, 4)];
					T1e = T1 - T2;
					T3 = T1 + T2;
					T1C = Ty - Tz;
					TA = Ty + Tz;
					T5 = Rm[WS(rs, 3)];
				   }
				   T1f = Ip[WS(rs, 4)];
				   T1g = Im[WS(rs, 3)];
				   {
					E Tb, Tc, TI, TJ;
					Tb = Rm[WS(rs, 1)];
					Tx = T4 - T5;
					T6 = T4 + T5;
					T1D = T1f - T1g;
					T1h = T1f + T1g;
					Tc = Rp[WS(rs, 6)];
					TI = Im[WS(rs, 1)];
					TJ = Ip[WS(rs, 6)];
					{
					     E T8, TH, TK, T9, TD, TE;
					     T8 = Rp[WS(rs, 2)];
					     Td = Tb + Tc;
					     TH = Tb - Tc;
					     T1A = TJ - TI;
					     TK = TI + TJ;
					     T9 = Rm[WS(rs, 5)];
					     TD = Ip[WS(rs, 2)];
					     TE = Im[WS(rs, 5)];
					     TL = TH + TK;
					     T1k = TH - TK;
					     Ta = T8 + T9;
					     TC = T8 - T9;
					     TF = TD + TE;
					     T1z = TD - TE;
					}
				   }
			      }
			      {
				   E T2E, TB, T1l, T1i, T3d, T3e, TM, T2F;
				   {
					E T7, TG, Te, T30, T31, T1j;
					T2E = T3 - T6;
					T7 = T3 + T6;
					T1j = TC - TF;
					TG = TC + TF;
					Te = Ta + Td;
					T2R = Ta - Td;
					TB = Tx + TA;
					T30 = TA - Tx;
					T31 = T1j - T1k;
					T1l = T1j + T1k;
					T1i = T1e - T1h;
					T3d = T1e + T1h;
					T20 = T7 - Te;
					Tf = T7 + Te;
					T3Q = FNMS(KP707106781, T31, T30);
					T32 = FMA(KP707106781, T31, T30);
					T3e = TG + TL;
					TM = TG - TL;
				   }
				   T3V = FMA(KP707106781, T3e, T3d);
				   T3f = FNMS(KP707106781, T3e, T3d);
				   T2a = FNMS(KP707106781, TM, TB);
				   TN = FMA(KP707106781, TM, TB);
				   T2F = T1A - T1z;
				   T1B = T1z + T1A;
				   T2f = FNMS(KP707106781, T1l, T1i);
				   T1m = FMA(KP707106781, T1l, T1i);
				   T3G = T2E - T2F;
				   T2G = T2E + T2F;
				   T2S = T1C - T1D;
				   T1E = T1C + T1D;
			      }
			 }
			 {
			      E T34, TS, T2H, Tm, T1u, T2I, T33, TX, Tq, T14, Tp, T1v, T12, Tr, T15;
			      E T16;
			      {
				   E Tj, TT, Ti, T1s, TR, Tk, TU, TV;
				   {
					E Tg, Th, TP, TQ;
					Tg = Rp[WS(rs, 1)];
					T3L = T2S - T2R;
					T2T = T2R + T2S;
					T26 = T1E - T1B;
					T1F = T1B + T1E;
					Th = Rm[WS(rs, 6)];
					TP = Ip[WS(rs, 1)];
					TQ = Im[WS(rs, 6)];
					Tj = Rp[WS(rs, 5)];
					TT = Tg - Th;
					Ti = Tg + Th;
					T1s = TP - TQ;
					TR = TP + TQ;
					Tk = Rm[WS(rs, 2)];
					TU = Ip[WS(rs, 5)];
					TV = Im[WS(rs, 2)];
				   }
				   {
					E Tn, To, T10, T11;
					Tn = Rm[0];
					{
					     E TO, Tl, T1t, TW;
					     TO = Tj - Tk;
					     Tl = Tj + Tk;
					     T1t = TU - TV;
					     TW = TU + TV;
					     T34 = TR - TO;
					     TS = TO + TR;
					     T2H = Ti - Tl;
					     Tm = Ti + Tl;
					     T1u = T1s + T1t;
					     T2I = T1s - T1t;
					     T33 = TT + TW;
					     TX = TT - TW;
					     To = Rp[WS(rs, 7)];
					}
					T10 = Im[0];
					T11 = Ip[WS(rs, 7)];
					Tq = Rp[WS(rs, 3)];
					T14 = Tn - To;
					Tp = Tn + To;
					T1v = T11 - T10;
					T12 = T10 + T11;
					Tr = Rm[WS(rs, 4)];
					T15 = Ip[WS(rs, 3)];
					T16 = Im[WS(rs, 4)];
				   }
			      }
			      {
				   E T13, T1x, T18, T35, T3g, T3h, T38, TY, T19;
				   {
					E T2U, T2J, T37, Tt, T36, T2V, T2M, T2K, T2L;
					T2U = T2H + T2I;
					T2J = T2H - T2I;
					{
					     E TZ, Ts, T1w, T17;
					     TZ = Tq - Tr;
					     Ts = Tq + Tr;
					     T1w = T15 - T16;
					     T17 = T15 + T16;
					     T37 = TZ + T12;
					     T13 = TZ - T12;
					     T2K = Tp - Ts;
					     Tt = Tp + Ts;
					     T1x = T1v + T1w;
					     T2L = T1v - T1w;
					     T36 = T14 + T17;
					     T18 = T14 - T17;
					}
					T2V = T2L - T2K;
					T2M = T2K + T2L;
					T3M = T2J - T2M;
					T2N = T2J + T2M;
					T3H = T2V - T2U;
					T2W = T2U + T2V;
					T35 = FMA(KP414213562, T34, T33);
					T3g = FNMS(KP414213562, T33, T34);
					T25 = Tm - Tt;
					Tu = Tm + Tt;
					T3h = FNMS(KP414213562, T36, T37);
					T38 = FMA(KP414213562, T37, T36);
				   }
				   T1n = FNMS(KP414213562, TS, TX);
				   TY = FMA(KP414213562, TX, TS);
				   T19 = FNMS(KP414213562, T18, T13);
				   T1o = FMA(KP414213562, T13, T18);
				   T3R = T3h - T3g;
				   T3i = T3g + T3h;
				   T2g = T19 - TY;
				   T1a = TY + T19;
				   T21 = T1x - T1u;
				   T1y = T1u + T1x;
				   T3W = T35 + T38;
				   T39 = T35 - T38;
			      }
			 }
		    }
		    {
			 E T27, T22, T2c, T2u, T2x, T2h, T2s, T2A, T2w, T2B, T2v;
			 {
			      E T1K, Tv, T1G, T1N, T1Q, T1b, T2b, T1p, Tw, T1d;
			      T1K = Tf - Tu;
			      Tv = Tf + Tu;
			      T1G = T1y + T1F;
			      T1N = T1F - T1y;
			      T1Q = FNMS(KP923879532, T1a, TN);
			      T1b = FMA(KP923879532, T1a, TN);
			      T2b = T1n - T1o;
			      T1p = T1n + T1o;
			      Tw = W[0];
			      T1d = W[1];
			      {
				   E T1T, T1O, T1W, T1S, T1X, T1R;
				   {
					E T1J, T1M, T1L, T1V, T1P, T1q;
					T1T = FNMS(KP923879532, T1p, T1m);
					T1q = FMA(KP923879532, T1p, T1m);
					{
					     E T1c, T1I, T1H, T1r;
					     T1c = Tw * T1b;
					     T1J = W[14];
					     T1H = Tw * T1q;
					     T1r = FMA(T1d, T1q, T1c);
					     T1M = W[15];
					     T1L = T1J * T1K;
					     T1I = FNMS(T1d, T1b, T1H);
					     Rm[0] = Tv + T1r;
					     Rp[0] = Tv - T1r;
					     T1V = T1M * T1K;
					     Im[0] = T1I - T1G;
					     Ip[0] = T1G + T1I;
					     T1P = W[16];
					}
					T1O = FNMS(T1M, T1N, T1L);
					T1W = FMA(T1J, T1N, T1V);
					T1S = W[17];
					T1X = T1P * T1T;
					T1R = T1P * T1Q;
				   }
				   {
					E T2r, T2n, T2q, T2p, T2z, T2t, T2o, T1Y, T1U;
					T27 = T25 + T26;
					T2r = T26 - T25;
					T2o = T20 - T21;
					T22 = T20 + T21;
					T1Y = FNMS(T1S, T1Q, T1X);
					T1U = FMA(T1S, T1T, T1R);
					T2n = W[22];
					T2q = W[23];
					Im[WS(rs, 4)] = T1Y - T1W;
					Ip[WS(rs, 4)] = T1W + T1Y;
					Rm[WS(rs, 4)] = T1O + T1U;
					Rp[WS(rs, 4)] = T1O - T1U;
					T2p = T2n * T2o;
					T2z = T2q * T2o;
					T2c = FMA(KP923879532, T2b, T2a);
					T2u = FNMS(KP923879532, T2b, T2a);
					T2x = FNMS(KP923879532, T2g, T2f);
					T2h = FMA(KP923879532, T2g, T2f);
					T2t = W[24];
					T2s = FNMS(T2q, T2r, T2p);
					T2A = FMA(T2n, T2r, T2z);
					T2w = W[25];
					T2B = T2t * T2x;
					T2v = T2t * T2u;
				   }
			      }
			 }
			 {
			      E T28, T2k, T2e, T2l, T2d;
			      {
				   E T1Z, T24, T23, T2j, T29, T2C, T2y;
				   T2C = FNMS(T2w, T2u, T2B);
				   T2y = FMA(T2w, T2x, T2v);
				   T1Z = W[6];
				   T24 = W[7];
				   Im[WS(rs, 6)] = T2C - T2A;
				   Ip[WS(rs, 6)] = T2A + T2C;
				   Rm[WS(rs, 6)] = T2s + T2y;
				   Rp[WS(rs, 6)] = T2s - T2y;
				   T23 = T1Z * T22;
				   T2j = T24 * T22;
				   T29 = W[8];
				   T28 = FNMS(T24, T27, T23);
				   T2k = FMA(T1Z, T27, T2j);
				   T2e = W[9];
				   T2l = T29 * T2h;
				   T2d = T29 * T2c;
			      }
			      {
				   E T4a, T4d, T3O, T4h, T4b, T40, T48, T4g, T4c, T42, T3Y;
				   {
					E T3N, T47, T43, T46, T3F, T45, T4f, T3K, T3J, T3S, T3X, T3Z, T49, T41, T3T;
					E T3U;
					{
					     E T44, T3I, T2m, T2i, T3P;
					     T44 = FNMS(KP707106781, T3H, T3G);
					     T3I = FMA(KP707106781, T3H, T3G);
					     T2m = FNMS(T2e, T2c, T2l);
					     T2i = FMA(T2e, T2h, T2d);
					     T3N = FMA(KP707106781, T3M, T3L);
					     T47 = FNMS(KP707106781, T3M, T3L);
					     Im[WS(rs, 2)] = T2m - T2k;
					     Ip[WS(rs, 2)] = T2k + T2m;
					     Rm[WS(rs, 2)] = T28 + T2i;
					     Rp[WS(rs, 2)] = T28 - T2i;
					     T43 = W[26];
					     T46 = W[27];
					     T3F = W[10];
					     T45 = T43 * T44;
					     T4f = T46 * T44;
					     T3K = W[11];
					     T3J = T3F * T3I;
					     T4a = FNMS(KP923879532, T3R, T3Q);
					     T3S = FMA(KP923879532, T3R, T3Q);
					     T3X = FNMS(KP923879532, T3W, T3V);
					     T4d = FMA(KP923879532, T3W, T3V);
					     T3Z = T3K * T3I;
					     T3P = W[12];
					     T49 = W[28];
					     T41 = T3P * T3X;
					     T3T = T3P * T3S;
					}
					T3O = FNMS(T3K, T3N, T3J);
					T4h = T49 * T4d;
					T4b = T49 * T4a;
					T40 = FMA(T3F, T3N, T3Z);
					T3U = W[13];
					T48 = FNMS(T46, T47, T45);
					T4g = FMA(T43, T47, T4f);
					T4c = W[29];
					T42 = FNMS(T3U, T3S, T41);
					T3Y = FMA(T3U, T3X, T3T);
				   }
				   {
					E T3t, T2X, T3p, T3s, T2D, T3r, T3B, T2Q, T2P, T3a, T3j, T3l, T3v, T3n, T3b;
					E T3c;
					{
					     E T2O, T3q, T4i, T4e, T2Z;
					     T4i = FNMS(T4c, T4a, T4h);
					     T4e = FMA(T4c, T4d, T4b);
					     Im[WS(rs, 3)] = T42 - T40;
					     Ip[WS(rs, 3)] = T40 + T42;
					     Rm[WS(rs, 3)] = T3O + T3Y;
					     Rp[WS(rs, 3)] = T3O - T3Y;
					     Im[WS(rs, 7)] = T4i - T4g;
					     Ip[WS(rs, 7)] = T4g + T4i;
					     Rm[WS(rs, 7)] = T48 + T4e;
					     Rp[WS(rs, 7)] = T48 - T4e;
					     T3t = FNMS(KP707106781, T2W, T2T);
					     T2X = FMA(KP707106781, T2W, T2T);
					     T2O = FMA(KP707106781, T2N, T2G);
					     T3q = FNMS(KP707106781, T2N, T2G);
					     T3p = W[18];
					     T3s = W[19];
					     T2D = W[2];
					     T3r = T3p * T3q;
					     T3B = T3s * T3q;
					     T2Q = W[3];
					     T2P = T2D * T2O;
					     T3a = FMA(KP923879532, T39, T32);
					     T3w = FNMS(KP923879532, T39, T32);
					     T3z = FMA(KP923879532, T3i, T3f);
					     T3j = FNMS(KP923879532, T3i, T3f);
					     T3l = T2Q * T2O;
					     T2Z = W[4];
					     T3v = W[20];
					     T3n = T2Z * T3j;
					     T3b = T2Z * T3a;
					}
					T2Y = FNMS(T2Q, T2X, T2P);
					T3D = T3v * T3z;
					T3x = T3v * T3w;
					T3m = FMA(T2D, T2X, T3l);
					T3c = W[5];
					T3u = FNMS(T3s, T3t, T3r);
					T3C = FMA(T3p, T3t, T3B);
					T3y = W[21];
					T3o = FNMS(T3c, T3a, T3n);
					T3k = FMA(T3c, T3j, T3b);
				   }
			      }
			 }
		    }
	       }
	       T3E = FNMS(T3y, T3w, T3D);
	       T3A = FMA(T3y, T3z, T3x);
	       Im[WS(rs, 1)] = T3o - T3m;
	       Ip[WS(rs, 1)] = T3m + T3o;
	       Rm[WS(rs, 1)] = T2Y + T3k;
	       Rp[WS(rs, 1)] = T2Y - T3k;
	       Im[WS(rs, 5)] = T3E - T3C;
	       Ip[WS(rs, 5)] = T3C + T3E;
	       Rm[WS(rs, 5)] = T3u + T3A;
	       Rp[WS(rs, 5)] = T3u - T3A;
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 16},
     {TW_NEXT, 1, 0}
};

static const hc2c_desc desc = { 16, "hc2cbdft2_16", twinstr, &GENUS, {136, 30, 70, 0} };

void X(codelet_hc2cbdft2_16) (planner *p) {
     X(khc2c_register) (p, hc2cbdft2_16, &desc, HC2C_VIA_DFT);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2cdft.native -compact -variables 4 -pipeline-latency 4 -sign 1 -n 16 -dif -name hc2cbdft2_16 -include hc2cb.h */

/*
 * This function contains 206 FP additions, 84 FP multiplications,
 * (or, 168 additions, 46 multiplications, 38 fused multiply/add),
 * 60 stack variables, 3 constants, and 64 memory accesses
 */
#include "hc2cb.h"

static void hc2cbdft2_16(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 30); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 30, MAKE_VOLATILE_STRIDE(64, rs)) {
	       E TB, T2L, T30, T1n, Tf, T1U, T2H, T3p, T1E, T1Z, TM, T31, T2s, T3k, T1i;
	       E T2M, Tu, T1Y, T2Q, T2X, T2T, T2Y, TY, T1d, T19, T1e, T2v, T2C, T2y, T2D;
	       E T1x, T1V;
	       {
		    E T3, T1j, TA, T1B, T6, Tx, T1m, T1C, Ta, TC, TF, T1y, Td, TH, TK;
		    E T1z;
		    {
			 E T1, T2, Ty, Tz;
			 T1 = Rp[0];
			 T2 = Rm[WS(rs, 7)];
			 T3 = T1 + T2;
			 T1j = T1 - T2;
			 Ty = Ip[0];
			 Tz = Im[WS(rs, 7)];
			 TA = Ty + Tz;
			 T1B = Ty - Tz;
		    }
		    {
			 E T4, T5, T1k, T1l;
			 T4 = Rp[WS(rs, 4)];
			 T5 = Rm[WS(rs, 3)];
			 T6 = T4 + T5;
			 Tx = T4 - T5;
			 T1k = Ip[WS(rs, 4)];
			 T1l = Im[WS(rs, 3)];
			 T1m = T1k + T1l;
			 T1C = T1k - T1l;
		    }
		    {
			 E T8, T9, TD, TE;
			 T8 = Rp[WS(rs, 2)];
			 T9 = Rm[WS(rs, 5)];
			 Ta = T8 + T9;
			 TC = T8 - T9;
			 TD = Ip[WS(rs, 2)];
			 TE = Im[WS(rs, 5)];
			 TF = TD + TE;
			 T1y = TD - TE;
		    }
		    {
			 E Tb, Tc, TI, TJ;
			 Tb = Rm[WS(rs, 1)];
			 Tc = Rp[WS(rs, 6)];
			 Td = Tb + Tc;
			 TH = Tb - Tc;
			 TI = Im[WS(rs, 1)];
			 TJ = Ip[WS(rs, 6)];
			 TK = TI + TJ;
			 T1z = TJ - TI;
		    }
		    {
			 E T7, Te, TG, TL;
			 TB = Tx + TA;
			 T2L = TA - Tx;
			 T30 = T1j + T1m;
			 T1n = T1j - T1m;
			 T7 = T3 + T6;
			 Te = Ta + Td;
			 Tf = T7 + Te;
			 T1U = T7 - Te;
			 {
			      E T2F, T2G, T1A, T1D;
			      T2F = Ta - Td;
			      T2G = T1B - T1C;
			      T2H = T2F + T2G;
			      T3p = T2G - T2F;
			      T1A = T1y + T1z;
			      T1D = T1B + T1C;
			      T1E = T1A + T1D;
			      T1Z = T1D - T1A;
			 }
			 TG = TC + TF;
			 TL = TH + TK;
			 TM = KP707106781 * (TG - TL);
			 T31 = KP707106781 * (TG + TL);
			 {
			      E T2q, T2r, T1g, T1h;
			      T2q = T3 - T6;
			      T2r = T1z - T1y;
			      T2s = T2q + T2r;
			      T3k = T2q - T2r;
			      T1g = TC - TF;
			      T1h = TH - TK;
			      T1i = KP707106781 * (T1g + T1h);
			      T2M = KP707106781 * (T1g - T1h);
			 }
		    }
	       }
	       {
		    E Ti, TT, TR, T1r, Tl, TO, TW, T1s, Tp, T14, T12, T1u, Ts, TZ, T17;
		    E T1v;
		    {
			 E Tg, Th, TP, TQ;
			 Tg = Rp[WS(rs, 1)];
			 Th = Rm[WS(rs, 6)];
			 Ti = Tg + Th;
			 TT = Tg - Th;
			 TP = Ip[WS(rs, 1)];
			 TQ = Im[WS(rs, 6)];
			 TR = TP + TQ;
			 T1r = TP - TQ;
		    }
		    {
			 E Tj, Tk, TU, TV;
			 Tj = Rp[WS(rs, 5)];
			 Tk = Rm[WS(rs, 2)];
			 Tl = Tj + Tk;
			 TO = Tj - Tk;
			 TU = Ip[WS(rs, 5)];
			 TV = Im[WS(rs, 2)];
			 TW = TU + TV;
			 T1s = TU - TV;
		    }
		    {
			 E Tn, To, T10, T11;
			 Tn = Rm[0];
			 To = Rp[WS(rs, 7)];
			 Tp = Tn + To;
			 T14 = Tn - To;
			 T10 = Im[0];
			 T11 = Ip[WS(rs, 7)];
			 T12 = T10 + T11;
			 T1u = T11 - T10;
		    }
		    {
			 E Tq, Tr, T15, T16;
			 Tq = Rp[WS(rs, 3)];
			 Tr = Rm[WS(rs, 4)];
			 Ts = Tq + Tr;
			 TZ = Tq - Tr;
			 T15 = Ip[WS(rs, 3)];
			 T16 = Im[WS(rs, 4)];
			 T17 = T15 + T16;
			 T1v = T15 - T16;
		    }
		    {
			 E Tm, Tt, T2O, T2P;
			 Tm = Ti + Tl;
			 Tt = Tp + Ts;
			 Tu = Tm + Tt;
			 T1Y = Tm - Tt;
			 T2O = TR - TO;
			 T2P = TT + TW;
			 T2Q = FMA(KP382683432, T2O, KP923879532 * T2P);
			 T2X = FNMS(KP923879532, T2O, KP382683432 * T2P);
		    }
		    {
			 E T2R, T2S, TS, TX;
			 T2R = TZ + T12;
			 T2S = T14 + T17;
			 T2T = FMA(KP382683432, T2R, KP923879532 * T2S);
			 T2Y = FNMS(KP923879532, T2R, KP382683432 * T2S);
			 TS = TO + TR;
			 TX = TT - TW;
			 TY = FMA(KP923879532, TS, KP382683432 * TX);
			 T1d = FNMS(KP382683432, TS, KP923879532 * TX);
		    }
		    {
			 E T13, T18, T2t, T2u;
			 T13 = TZ - T12;
			 T18 = T14 - T17;
			 T19 = FNMS(KP382683432, T18, KP923879532 * T13);
			 T1e = FMA(KP382683432, T13, KP923879532 * T18);
			 T2t = Ti - Tl;
			 T2u = T1r - T1s;
			 T2v = T2t - T2u;
			 T2C = T2t + T2u;
		    }
		    {
			 E T2w, T2x, T1t, T1w;
			 T2w = Tp - Ts;
			 T2x = T1u - T1v;
			 T2y = T2w + T2x;
			 T2D = T2x - T2w;
			 T1t = T1r + T1s;
			 T1w = T1u + T1v;
			 T1x = T1t + T1w;
			 T1V = T1w - T1t;
		    }
	       }
	       {
		    E Tv, T1F, T1b, T1N, T1p, T1P, T1L, T1R;
		    Tv = Tf + Tu;
		    T1F = T1x + T1E;
		    {
			 E TN, T1a, T1f, T1o;
			 TN = TB + TM;
			 T1a = TY + T19;
			 T1b = TN + T1a;
			 T1N = TN - T1a;
			 T1f = T1d + T1e;
			 T1o = T1i + T1n;
			 T1p = T1f + T1o;
			 T1P = T1o - T1f;
			 {
			      E T1I, T1K, T1H, T1J;
			      T1I = Tf - Tu;
			      T1K = T1E - T1x;
			      T1H = W[14];
			      T1J = W[15];
			      T1L = FNMS(T1J, T1K, T1H * T1I);
			      T1R = FMA(T1J, T1I, T1H * T1K);
			 }
		    }
		    {
			 E T1q, T1G, Tw, T1c;
			 Tw = W[0];
			 T1c = W[1];
			 T1q = FMA(Tw, T1b, T1c * T1p);
			 T1G = FNMS(T1c, T1b, Tw * T1p);
			 Rp[0] = Tv - T1q;
			 Ip[0] = T1F + T1G;
			 Rm[0] = Tv + T1q;
			 Im[0] = T1G - T1F;
		    }
		    {
			 E T1Q, T1S, T1M, T1O;
			 T1M = W[16];
			 T1O = W[17];
			 T1Q = FMA(T1M, T1N, T1O * T1P);
			 T1S = FNMS(T1O, T1N, T1M * T1P);
			 Rp[WS(rs, 4)] = T1L - T1Q;
			 Ip[WS(rs, 4)] = T1R + T1S;
			 Rm[WS(rs, 4)] = T1L + T1Q;
			 Im[WS(rs, 4)] = T1S - T1R;
		    }
	       }
	       {
		    E T25, T2j, T29, T2l, T21, T2b, T2h, T2n;
		    {
			 E T23, T24, T27, T28;
			 T23 = TB - TM;
			 T24 = T1d - T1e;
			 T25 = T23 + T24;
			 T2j = T23 - T24;
			 T27 = T19 - TY;
			 T28 = T1n - T1i;
			 T29 = T27 + T28;
			 T2l = T28 - T27;
		    }
		    {
			 E T1W, T20, T1T, T1X;
			 T1W = T1U + T1V;
			 T20 = T1Y + T1Z;
			 T1T = W[6];
			 T1X = W[7];
			 T21 = FNMS(T1X, T20, T1T * T1W);
			 T2b = FMA(T1X, T1W, T1T * T20);
		    }
		    {
			 E T2e, T2g, T2d, T2f;
			 T2e = T1U - T1V;
			 T2g = T1Z - T1Y;
			 T2d = W[22];
			 T2f = W[23];
			 T2h = FNMS(T2f, T2g, T2d * T2e);
			 T2n = FMA(T2f, T2e, T2d * T2g);
		    }
		    {
			 E T2a, T2c, T22, T26;
			 T22 = W[8];
			 T26 = W[9];
			 T2a = FMA(T22, T25, T26 * T29);
			 T2c = FNMS(T26, T25, T22 * T29);
			 Rp[WS(rs, 2)] = T21 - T2a;
			 Ip[WS(rs, 2)] = T2b + T2c;
			 Rm[WS(rs, 2)] = T21 + T2a;
			 Im[WS(rs, 2)] = T2c - T2b;
		    }
		    {
			 E T2m, T2o, T2i, T2k;
			 T2i = W[24];
			 T2k = W[25];
			 T2m = FMA(T2i, T2j, T2k * T2l);
			 T2o = FNMS(T2k, T2j, T2i * T2l);
			 Rp[WS(rs, 6)] = T2h - T2m;
			 Ip[WS(rs, 6)] = T2n + T2o;
			 Rm[WS(rs, 6)] = T2h + T2m;
			 Im[WS(rs, 6)] = T2o - T2n;
		    }
	       }
	       {
		    E T2A, T38, T2I, T3a, T2V, T3d, T33, T3f, T2z, T2E;
		    T2z = KP707106781 * (T2v + T2y);
		    T2A = T2s + T2z;
		    T38 = T2s - T2z;
		    T2E = KP707106781 * (T2C + T2D);
		    T2I = T2E + T2H;
		    T3a = T2H - T2E;
		    {
			 E T2N, T2U, T2Z, T32;
			 T2N = T2L + T2M;
			 T2U = T2Q - T2T;
			 T2V = T2N + T2U;
			 T3d = T2N - T2U;
			 T2Z = T2X + T2Y;
			 T32 = T30 - T31;
			 T33 = T2Z + T32;
			 T3f = T32 - T2Z;
		    }
		    {
			 E T2J, T35, T34, T36;
			 {
			      E T2p, T2B, T2K, T2W;
			      T2p = W[2];
			      T2B = W[3];
			      T2J = FNMS(T2B, T2I, T2p * T2A);
			      T35 = FMA(T2B, T2A, T2p * T2I);
			      T2K = W[4];
			      T2W = W[5];
			      T34 = FMA(T2K, T2V, T2W * T33);
			      T36 = FNMS(T2W, T2V, T2K * T33);
			 }
			 Rp[WS(rs, 1)] = T2J - T34;
			 Ip[WS(rs, 1)] = T35 + T36;
			 Rm[WS(rs, 1)] = T2J + T34;
			 Im[WS(rs, 1)] = T36 - T35;
		    }
		    {
			 E T3b, T3h, T3g, T3i;
			 {
			      E T37, T39, T3c, T3e;
			      T37 = W[18];
			      T39 = W[19];
			      T3b = FNMS(T39, T3a, T37 * T38);
			      T3h = FMA(T39, T38, T37 * T3a);
			      T3c = W[20];
			      T3e = W[21];
			      T3g = FMA(T3c, T3d, T3e * T3f);
			      T3i = FNMS(T3e, T3d, T3c * T3f);
			 }
			 Rp[WS(rs, 5)] = T3b - T3g;
			 Ip[WS(rs, 5)] = T3h + T3i;
			 Rm[WS(rs, 5)] = T3b + T3g;
			 Im[WS(rs, 5)] = T3i - T3h;
		    }
	       }
	       {
		    E T3m, T3E, T3q, T3G, T3v, T3J, T3z, T3L, T3l, T3o;
		    T3l = KP707106781 * (T2D - T2C);
		    T3m = T3k + T3l;
		    T3E = T3k - T3l;
		    T3o = KP707106781 * (T2v - T2y);
		    T3q = T3o + T3p;
		    T3G = T3p - T3o;
		    {
			 E T3t, T3u, T3x, T3y;
			 T3t = T2L - T2M;
			 T3u = T2X - T2Y;
			 T3v = T3t + T3u;
			 T3J = T3t - T3u;
			 T3x = T31 + T30;
			 T3y = T2Q + T2T;
			 T3z = T3x - T3y;
			 T3L = T3y + T3x;
		    }
		    {
			 E T3r, T3B, T3A, T3C;
			 {
			      E T3j, T3n, T3s, T3w;
			      T3j = W[10];
			      T3n = W[11];
			      T3r = FNMS(T3n, T3q, T3j * T3m);
			      T3B = FMA(T3n, T3m, T3j * T3q);
			      T3s = W[12];
			      T3w = W[13];
			      T3A = FMA(T3s, T3v, T3w * T3z);
			      T3C = FNMS(T3w, T3v, T3s * T3z);
			 }
			 Rp[WS(rs, 3)] = T3r - T3A;
			 Ip[WS(rs, 3)] = T3B + T3C;
			 Rm[WS(rs, 3)] = T3r + T3A;
			 Im[WS(rs, 3)] = T3C - T3B;
		    }
		    {
			 E T3H, T3N, T3M, T3O;
			 {
			      E T3D, T3F, T3I, T3K;
			      T3D = W[26];
			      T3F = W[27];
			      T3H = FNMS(T3F, T3G, T3D * T3E);
			      T3N = FMA(T3F, T3E, T3D * T3G);
			      T3I = W[28];
			      T3K = W[29];
			      T3M = FMA(T3I, T3J, T3K * T3L);
			      T3O = FNMS(T3K, T3J, T3I * T3L);
			 }
			 Rp[WS(rs, 7)] = T3H - T3M;
			 Ip[WS(rs, 7)] = T3N + T3O;
			 Rm[WS(rs, 7)] = T3H + T3M;
			 Im[WS(rs, 7)] = T3O - T3N;
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 16},
     {TW_NEXT, 1, 0}
};

static const hc2c_desc desc = { 16, "hc2cbdft2_16", twinstr, &GENUS, {168, 46, 38, 0} };

void X(codelet_hc2cbdft2_16) (planner *p) {
     X(khc2c_register) (p, hc2cbdft2_16, &desc, HC2C_VIA_DFT);
}
#endif				/* HAVE_FMA */
