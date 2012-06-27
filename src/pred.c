#include <math.h>
void pred(double * Ehat, double * ghat, double * R, double * BFt,  double * GWG, double * y, int * nin, int * din, int * nxin)
{
  int n=*nin;
  int nx=*nxin;
  int d=*din;
  int i,j, ii, jj, id, jd, iid;
  double tot, gt;
  for(i=0; i < nx; i++)
  {
    gt=0;
    id=i*d;
    for(j=0; j < n; j++)
    {
      jd=j*d;
      tot=0;
      for(ii=0; ii<d;ii++)
      {
        iid=ii*d;
        for(jj=0; jj<d;jj++)
        {
          tot+= GWG[iid+jj] * (R[id + ii] - BFt[jd + ii])*(R[id+jj] - BFt[jd+jj]);
        }
      }
      ghat[j]=exp(-0.5*tot);     
      gt+=ghat[j];      
    }
    
    tot=0;
    for(j=0; j < n; j++)
    {
      tot+=y[j]*ghat[j]/gt; 
    }
    Ehat[i]=tot;
  }
}


			 
