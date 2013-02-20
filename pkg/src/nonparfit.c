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
