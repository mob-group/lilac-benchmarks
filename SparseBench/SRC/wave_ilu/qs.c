#include <stdlib.h>
#include <math.h>
#include "qs.h"

int qsi(int *jdata,int len) {return qsip(jdata,0,len);}

int qsip(int *jdata,int *jperm,int len)
{int l=0; return qsip1(jdata,jperm,&len,&l);}

int qsip1(int *jdata,int *jperm,int *len,int *l)
{
  int r = 0, w = 0, b = *len-1, mid = jdata[(int)(*len/2)],
    i,jt,pt,ll=*l+1,nl;

  /*  printf("Input:");for (i=0;i<*len;i++) printf("(%d,%e) ",jdata[i],rdata[i]); printf("\n");*/
  if (*len<=1) {
    return 0;
  } else if (*len==2) {
    if (jdata[1]<jdata[0]) {
      jt = jdata[1]; jdata[1] = jdata[0]; jdata[0] = jt;
      if (jperm) {
	pt = jperm[1]; jperm[1] = jperm[0]; jperm[0] = pt;}
    }
    return 0;
  } else {
    for (i=1; i<*len; i++) {if (jdata[i]<jdata[i-1]) goto sort;}
    return 0;
  }
 sort: /*printf("midvalue: %d\n",mid);*/
  for (;w<=b;) {
    if (jdata[w]>mid) {
      jt = jdata[w]; jdata[w] = jdata[b]; jdata[b] = jt;
      if (jperm) {
	pt = jperm[w]; jperm[w] = jperm[b]; jperm[b] = pt;}
      b--;
    }
    if (jdata[w]<mid) {
      if (r==w) {
	r++; w++;
      } else {
	jt = jdata[w]; jdata[w] = jdata[r]; jdata[r] = jt;
	if (jperm) {
	  pt = jperm[w]; jperm[w] = jperm[r]; jperm[r] = pt;}
	r++;
      }
    }
    if (jdata[w]==mid) {w++;}
  }
 done:
  /*  printf("Through:");for (i=0;i<*len;i++) printf("(%d,%e) ",jdata[i],rdata[i]); printf("\n");*/
  nl = r;        qsip1(jdata,jperm,&nl,&ll);
  nl = *len-b-1;
  {int *jpl=0; if (jperm) jpl=jperm+b+1; qsip1(jdata+b+1,jpl,&nl,&ll);}

  /*
  printf("Result%d: ",*l);
  for (i=0; i<*len; i++) printf("(%d,%e) ",jdata[i],rdata[i]);
  printf("\n");
  */
  return 0;
}

int qsir(int *jdata,double *rdata,int len)
{int l=0; return qsir1(jdata,rdata,&len,&l);}

int qsir1(int *jdata,double *rdata,int *len,int *l)
{
  double rt;
  int r = 0, w = 0, b = *len-1, mid = jdata[(int)(*len/2)],
    i,jt,ll=*l+1,nl;

  /*  printf("Input:");for (i=0;i<*len;i++) printf("(%d,%e) ",jdata[i],rdata[i]); printf("\n");*/
  if (*len<=1) {
    return 0;
  } else if (*len==2) {
    if (jdata[1]<jdata[0]) {
      jt = jdata[1]; jdata[1] = jdata[0]; jdata[0] = jt;
      rt = rdata[1]; rdata[1] = rdata[0]; rdata[0] = rt;
    }
    return 0;
  } else {
    for (i=1; i<*len; i++) {if (jdata[i]<jdata[i-1]) goto sort;}
    return 0;
  }
 sort: /*printf("midvalue: %d\n",mid);*/
  for (;w<=b;) {
    if (jdata[w]>mid) {
      jt = jdata[w]; jdata[w] = jdata[b]; jdata[b] = jt;
      rt = rdata[w]; rdata[w] = rdata[b]; rdata[b] = rt;
      b--;
    }
    if (jdata[w]<mid) {
      if (r==w) {
	r++; w++;
      } else {
	jt = jdata[w]; jdata[w] = jdata[r]; jdata[r] = jt;
	rt = rdata[w]; rdata[w] = rdata[r]; rdata[r] = rt;
	r++;
      }
    }
    if (jdata[w]==mid) {w++;}
  }
 done:
  /*  printf("Through:");for (i=0;i<*len;i++) printf("(%d,%e) ",jdata[i],rdata[i]); printf("\n");*/
  nl = r;        qsir1(jdata,rdata,&nl,&ll);
  nl = *len-b-1; qsir1(jdata+b+1,rdata+b+1,&nl,&ll);

  /*
  printf("Result%d: ",*l);
  for (i=0; i<*len; i++) printf("(%d,%e) ",jdata[i],rdata[i]);
  printf("\n");
  */
  return 0;
}

int qsijr(int *idata,int *jdata,double *rdata,int len)
{int l=0; return qsijr1(idata,jdata,rdata,&len,&l);}

int qsijr1(int *idata,int *jdata,double *rdata,int *len,int *l)
{
  double rt;
  int r = 0, w = 0, b = *len-1, mid = jdata[(int)(*len/2)],
    i,it,jt,ll=*l+1,nl;

  /*  printf("Input:");for (i=0;i<*len;i++) printf("(%d,%e) ",jdata[i],rdata[i]); printf("\n");*/
  if (*len<=1) {
    return 0;
  } else if (*len==2) {
    if (jdata[1]<jdata[0]) {
      it = idata[1]; idata[1] = idata[0]; idata[0] = it;
      jt = jdata[1]; jdata[1] = jdata[0]; jdata[0] = jt;
      rt = rdata[1]; rdata[1] = rdata[0]; rdata[0] = rt;
    }
    return 0;
  } else {
    for (i=1; i<*len; i++) {if (jdata[i]<jdata[i-1]) goto sort;}
    return 0;
  }
 sort: /*printf("midvalue: %d\n",mid);*/
  for (;w<=b;) {
    if (jdata[w]>mid) {
      it = idata[w]; idata[w] = idata[b]; idata[b] = it;
      jt = jdata[w]; jdata[w] = jdata[b]; jdata[b] = jt;
      rt = rdata[w]; rdata[w] = rdata[b]; rdata[b] = rt;
      b--;
    }
    if (jdata[w]<mid) {
      if (r==w) {
	r++; w++;
      } else {
	it = idata[w]; idata[w] = idata[r]; idata[r] = it;
	jt = jdata[w]; jdata[w] = jdata[r]; jdata[r] = jt;
	rt = rdata[w]; rdata[w] = rdata[r]; rdata[r] = rt;
	r++;
      }
    }
    if (jdata[w]==mid) {w++;}
  }
 done:
  /*  printf("Through:");for (i=0;i<*len;i++) printf("(%d,%e) ",jdata[i],rdata[i]); printf("\n");*/
  nl = r;        qsijr1(idata,jdata,rdata,&nl,&ll);
  nl = *len-b-1; qsijr1(idata+b+1,jdata+b+1,rdata+b+1,&nl,&ll);

  /*
  printf("Result%d: ",*l);
  for (i=0; i<*len; i++) printf("(%d,%e) ",jdata[i],rdata[i]);
  printf("\n");
  */
  return 0;
}

int qsi_(int *jdata,int *len)
{int l=0; return qsip1(jdata,0,len,&l);}
int QSI(int *jdata,int *len)
{int l=0; return qsip1(jdata,0,len,&l);}

int qsip_(int *jdata,int *jperm,int *len)
{int l=0; return qsip1(jdata,jperm,len,&l);}
int QSIP(int *jdata,int *jperm,int *len)
{int l=0; return qsip1(jdata,jperm,len,&l);}

int qsir_(int *jdata,double *rdata,int *len)
{int l=0; return qsir1(jdata,rdata,len,&l);}
int QSIR(int *jdata,double *rdata,int *len)
{int l=0; return qsir1(jdata,rdata,len,&l);}

int qsijr_(int *idata,int *jdata,double *rdata,int *len)
{int l=0; return qsijr1(idata,jdata,rdata,len,&l);}
int QSIJR(int *idata,int *jdata,double *rdata,int *len)
{int l=0; return qsijr1(idata,jdata,rdata,len,&l);}

