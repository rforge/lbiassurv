/* Author: Pierre-Jerome Bergeron 

This is the Vardi EM algorithm from Biometrika (1989).
Provided as is. If you publish something using this code, please make sure to mention me....
*/

//This little function counts the number of time each unique entry appears
//in a sorted vector

void countsorted (double *x, int *nx, double *y) {
  int i, j=0;
  
  y[0]=1;
  for (i=1; i<*nx; i++) {
    if (x[i]==x[i-1]) y[j]++;
    else {j++; y[j]=1;}
  }
}

//These functions are going to be used at some point

void insert_value(double *value, double *superset, int *nx, double *subset, int *ny, double *f, double *results) {
  int i, j=0;
  
  for (i=0; i<*nx; i++) {
    if (j<*ny) {
      if (superset[i]==subset[j]) {
        results[i]=f[j];
        j++;
      }
      else results[i]=value[0];
    }
    else results[i]=value[0];
  }
}

//This is the Vardi algorithm coded in C
void vardir(double *x, double *y, int *m, int *n, double *tj, int *h, 
            int *maxiter, double *tol, double *pvec, int *iter, 
            double *xu, int *mx, double *yu, int *ny){
  int i,j, finished=0;
  double pold[h[0]],pnew[*h],val=0.0;
  double xcounts[mx[0]], ycounts[ny[0]],txcounts[h[0]], tycounts[h[0]];
  double tcounts[h[0]], tpk[h[0]],temp,tmp;
  
  for (i=0;i<*mx;i++) {
    xcounts[i]=0.0;
  }
  for (i=0;i<*ny;i++) {
    ycounts[i]=0.0;
  }
  for (i=0;i<*h;i++) {
    txcounts[i]=0.0;
    tycounts[i]=0.0;
  }
  /*Initializing xcounts, ycounts, tcounts, pold */
  countsorted(x,m,xcounts);
  countsorted(y,n,ycounts);
  iter[0]=0;    
  insert_value(&val,tj,h,xu,mx,xcounts,txcounts);
  insert_value(&val,tj,h,yu,ny,ycounts,tycounts);
  for (i=0;i<*h;i++) {
    tcounts[i]=txcounts[i]+tycounts[i];
    pold[i]=tcounts[i]/h[0];
  }
  
  while ((iter[0]<maxiter[0]) && (!finished)) {
    temp=0.0;
    for (j=h[0];j>0;j--) {
      temp+=pold[j-1]/tj[j-1];
      tpk[j-1]=tycounts[j-1]/temp;
    } //for j in h:1
    temp=0.0;
    for (j=0;j<*h;j++) {
      temp+=tpk[j];
      pnew[j]=(1.0/(m[0]+n[0]))*(txcounts[j]+(1.0/tj[j])*pold[j]*temp);
    }         

    tmp=0.0;
    for (j=0;j<*h;j++) {
      if (pold[j]<pnew[j])
        tmp+=pnew[j]-pold[j];
      else
        tmp+=pold[j]-pnew[j];
    }
    if (tmp<*tol) finished=1;
    for (j=0;j<*h;j++) {
      pold[j]=pnew[j];
    }
    iter[0]++;
  } //while
  for (j=0;j<h[0];j++) pvec[j]=pnew[j];
}



//This is a clumsily done homemade KM survival estimator. 
//Appears to work properly
void kaplanmeierr(double *t, int *d, int *n, double *ut, int *un, 
                  int *di, int *ni, double *St) {
  int i=0,j=0,k, ri=*n, fi;
  double curtime=ut[0], cursurv=1.0;

  while (t[i]<curtime) {
    i++;
    ri--;
  }
  for  (j=0;j<*un;j++) {
    curtime=ut[j];
    k=0;
    fi=0;
    while ((t[i]==curtime) && (i<*n)) {
      if (d[i]==1) fi++;
      k++;
      i++;
    }
    di[j]=fi;
    ni[j]=ri;
    cursurv *= (1.0-((double) fi/ (double) ri));
    St[j]=cursurv;
    ri-=k;  
    if (j<=(*un-1))
      while ((i<*n)&&(t[i]!=ut[j+1])) {
        i++;
  ri--;
      }
  }
}


void trunckaplanmeierr(double *ti, double *te, int *d, int *n, double *ute, int *une, int *di, int *ni, double *St) {
  int i1=0, i2=0,j=0,k, ri=0, fi;
  double curtime=ute[0], cursurv=1.0;

  while ((i1<*n)&&(ti[i1]<te[i2])) {
    ri++;
    i1++;
  }
  while (te[i2]<curtime) {
    i2++;
    ri--;  
    while ((i1<*n)&&(ti[i1]<te[i2])) {
      i1++;
      ri++;
    }
    
  }
  for  (j=0;j<*une;j++) {
    curtime=ute[j];
    k=0;
    fi=0;
    while ((te[i2]==curtime) && (i2<*n)) {
      if (d[i2]==1) fi++;
      k++;
      i2++;
    }
    di[j]=fi;
    ni[j]=ri;
    if (cursurv>0)
      cursurv *= (1.0-((double) fi/ (double) ri));
    St[j]=cursurv;
    ri-=k;  
    if (j<=(*une-1))
      while ((i1<*n)&&(ti[i1]<te[i2])) {
        i1++;
        ri++;
      }
      while ((i2<*n)&&(te[i2]!=ute[j+1])) {
        i2++;
	ri--;
        while ((i1<*n)&&(ti[i1]<te[i2])) {
          i1++;
          ri++;
        }
      }
  }
}






void weiwnr(double *x, double *y, double *d, int *n, double *p, double *wn, double *s0sq) {
  int i=0, j=0, k=0;
  double temp=0.0, n2=((double) n[0])*((double) n[0]), n3=n2*((double) n[0]);
  double sum1, sum2, sum3, sum5, sum8;
  
  for (i=0;i<n[0];i++) {
    for (j=0; j<n[0]; j++) {
      temp += ((double) ((x[i]>y[j])*d[j])) - ((double) (x[i]<y[j]))-p[0];
    }   
  }
  wn[0] = temp/n2;
  temp=0.0;
  for (j=0;j<n[0];j++) {
    sum1=0.0;
    sum2=0.0;
    sum3=0.0;
    sum5=0.0;
    sum8=0.0;   
    for (i=0;i<n[0];i++) {
      sum1+=(double) (y[i]>x[j]);
      sum2+=(double) (x[i]>y[j]);
      sum3+=(double) ((y[i]<=x[j])*d[i]);
      sum5+=(double) (x[i]<=y[j]);
      sum8+=(double) ((y[i]<=x[j])*d[i]);
    }
    temp+=sum1*sum1+sum2*sum2*d[j] + 2.0*sum3*sum2*d[j]+2.0*sum5*sum1 +2.0*sum5*sum8 -2.0*sum2*sum1*d[j];    
  }
  s0sq[0]=temp/n3;
}


void justweiwnr(double *x, double *y, double *d, int *n, double *p, double *wn) {
  int i=0, j=0, k=0;
  double temp=0.0, n2=((double) n[0])*((double) n[0]);
  
  
  for (i=0;i<n[0];i++) {
    for (j=0; j<n[0]; j++) {
      temp += ((double) ((x[i]>y[j])*d[j])) - ((double) (x[i]<y[j]))-p[0];
    }   
  }
  wn[0] = temp/n2;
}
