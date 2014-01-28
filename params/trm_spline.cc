template <class X> void splint(Subs::Buffer1D<X>& xa, Subs::Buffer1D<X>& ya, Subs::Buffer1D<X>& y2a, const X& x, X &y)
{
	int k;
        X h,b,a;

	int n=xa.size();
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) 
	  throw Subs::Subs_Error("splint(Subs::Buffer1D<X>& xa, Subs::Buffer1D<X>& ya, Subs::Buffer1D<X>& y2a, const X& x, X &y): Bad xa input");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]
		+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

template <class X> void spline(Subs::Buffer1D<X>& x, Subs::Buffer1D<X>& y, const X& yp1, const X& ypn,
	Buffer1D<X>& y2)
{
	int i,k;
	X p,qn,sig,un;

	int n=y2.size();
	Subs::Buffer1D<X> u(n-1);
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}

template <class X> void spline(const Subs::Array1D<X>&  x, const Subs::Array1D<X>&  y,
	    unsigned long size, Subs::Array1D<X>& xint, 
	    Subs::Array1D<X>& yint){

  if(x.size() != y.size()) 
    throw Subs::Subs_Error("Error returned inside spline: x and y arrays differ in size\n");
  Subs::Buffer1D<X> vecx(x.size()), vecy(x.size()), d2y(x.size());
  for(int i=0; i< x.size(); i++){vecx[i] = x[i]; vecy[i]=y[i];}
  X yp1, ypn;
  yp1 = 1.0e30; ypn = 1.0e30;
  spline(vecx,vecy,yp1,ypn,d2y);
  for(int i=0; i<size; i++){
    X xtmp, ytmp;
    xtmp = x.min() + (x.max()-x.min())*i/size;
    splint(vecx,vecy,d2y,xtmp,ytmp);
    xint.push_back(xtmp);
    yint.push_back(ytmp);
  }      
}
