#include <stdio.h>
#include <math.h>
#include "path.h"

double ran3()
{
  static int ma[60],mj,mk,mz,i,ii,k,inext,inextp,iff,mbig,mseed;
  
  
  if (iff ==0) {
    iff=1;
    mbig = 1000000000;
    mseed = 161803398;
    mz =0;
    mj = mseed ;
    mj = mj % mbig;
    ma[55] = mj;
    mk = 1;
    for (i=1;i<55;i++){
      ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk<mz) mk=mk+mbig;
      mj = ma[ii];
    }
    for(k=1;k<=4;k++){
      for(i=1;i<=55;i++){
	ma[i] = ma[i] - ma[1 + ((i+30)%55)];
	if (ma[i]<mz) ma[i] = ma[i] + mbig;
      }
    }
    inext=0;
    inextp=31;
   }
  if (++inext == 56)   inext  =1;
  if (++inextp ==56)  inextp =1;
  mj = ma[inext] - ma[inextp];
  if (mj<mz) mj=mj+mbig;
  ma[inext]=mj;
  return  (double)mj / mbig;
}
    


double gssran()

{
  double fac,rsq,v1,v2;
  static int iset=0;
  static double gset;

  if (iset==0)  {
    do{
      v1 =2.*ran3() -1.;
      v2 =2.*ran3() -1.;
      rsq = v1*v1 + v2*v2;
    } while ((rsq >= 1)|| (rsq ==0));
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset =1;
    return v2*fac;
  } else {
    iset =0;
    return gset;
  }
}




void gssbivar (double *r,
               double *v)
{
    
    double fac,rsq,v1,v2;

    /* draw bivariate random numbers */

    do{
      v1 =2.*ran3() -1.;
      v2 =2.*ran3() -1.;
      rsq = v1*v1 + v2*v2;
    } while ((rsq >= 1)|| (rsq ==0));
    fac = sqrt(-2.*log(rsq)/rsq);
   
    *r = v1*fac * bivar.sqrts11;
    *v = v2*fac * bivar.sqrtSos11 + (*r) * bivar.s12os11;
}



void quicksort(double *x, int first, int last)
{
  int pivot,j,i;
  double temp;
  
  if(first<last){
    pivot=first;
    i=first;
    j=last;

    while(i<j){
      while(x[i]<=x[pivot]&&i<last)
	i++;
      while(x[j]>x[pivot])
	j--;
      if(i<j){
	temp=x[i];
	x[i]=x[j];
	x[j]=temp;
      }
    }
    
    temp=x[pivot];
    x[pivot]=x[j];
    x[j]=temp;
    quicksort(x,first,j-1);
    quicksort(x,j+1,last);
    
  }
}
