#Author: Pierre-Jerome Bergeron 
# 
# Based on Vardi, Y. (1989) Multiplicative censoring, etc... Biometrika
# Provided as is. If you publish something using this code, make sure to acknowledge me...

#dyn.load("Dropbox/pjerome-vahid/ccode/pjbCfun.so")

#This is the Vardi algorithm, it takes data of the form (T, delta)
#Appears pretty robust
vardifit<-function(tt,d,maxiter=500,tol=1e-6) {
  
  #Counts the number of times each unique element appears in a sorted vector
  svcounts<-function(x) {
    .C("countsorted",PACKAGE="lbiassurv",as.double(x),as.integer(length(x)),y=double(length(unique(x))))$y
  }
  
  #Removes zeros from a vector
  remzero<-function(x) {x[x!=0]}
  
  #Function to split time and failure indicators into two seperate vectors of
  #uncensored times and censored times, the way Vardi wrote his algorithm
  splitxy<-function(tt,delta) {
    x<-tt[delta==1]
    y<-tt[delta==0]
    list(x=x,y=y) 
  }
  
# The algorithm starts from here
  ot<-order(tt)
  xy<-splitxy(tt[ot],d[ot])
  tj<-unique(tt[ot])
  h<-length(tj)
  m<-length(xy$x)
  n<-length(xy$y)
  iter=0
  if (n==0) {
    if (h==m) 
      vout<-list(pvec=rep(1/h,h),iter=iter)
    else
      vout<-list(pvec=svcounts(xy$x)/m,iter=iter)
  }
  else if (h==(m+n)) {
    vout<-.C("vardir",PACKAGE="lbiassurv",as.double(xy$x),as.double(xy$y),as.integer(m),as.integer(n),
             as.double(tj),as.integer(h),as.integer(maxiter), as.double(tol),
             pvec=double(h),iter=integer(1),as.double(xy$x),
             as.integer(m),as.double(xy$y),as.integer(n))
  }
  else { 
    xu<-unique(xy$x)
    yu<-unique(xy$y)
    vout<-.C("vardir",PACKAGE="lbiassurv",as.double(xy$x),as.double(xy$y),as.integer(m),as.integer(n),
             as.double(tj),as.integer(h),as.integer(maxiter), as.double(tol),
             pvec=double(h),iter=integer(1),as.double(xu),
             as.integer(length(xu)),as.double(yu),as.integer(length(yu)))
  }
  ret<-list(times=tj,lbpvec=vout$pvec,ubpvec=(vout$pvec/tj)/sum(vout$pvec/tj),
            iterations=vout$iter, conv=(iter<maxiter))
  return(ret)
}

# This function generates data with length-baised censoring 

lbiasgenerate<-function(size)
{
  lbtimevec=rgamma(size,shape=2,scale=1)
  lbtruncvec=runif(size,max=lbtimevec)
  residtimevec=lbtimevec-lbtruncvec
  residcensvec=rexp(size,rate=.5)
  obstimevec=lbtruncvec+pmin(residtimevec,residcensvec)
  deltavec=1*(residtimevec<residcensvec)
  return(list(data=obstimevec,censor=deltavec))
}


