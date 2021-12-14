#include <iostream> 
#include <vector> 
#include <fstream> 
#include <time.h> 
#include <Cantera.h> 
#include <onedim.h> 
#include <IdealGasMix.h> 
#include <equilibrium.h> 
#include <transport.h> 
using namespace std; 
using namespace Cantera; 
double Combust(vector<double> vars); 
int RSM(); 
int RSM2D(); 
double GRGM(vector<double> &x0, double f0, double alpha, double lcon, double rcon, int dim); 
double norm(vector<double> x); 
double inprod(vector<double>x, vector<double>y); 
vector<double> Mv(double **M,vector<double> v,double row, double col); 
double **MM(double **M1, double **M2,double a, double b, double c); 
vector<double> LUSolve(double **M,vector<double> v); 
double sum(vector<double> x); 
double newfun(vector<double> x, vector<double> b, double &g, double &H); 
double newfun2D(vector<double> x, vector<double> b, vector<double> &grad, double **&Hess); 


double fun2(double x); 
double fun2D(double x,double y); 
double **MT(double **M,double n, double p); 
double max(double x, double y); 
double min(double x, double y); 
int sign(double x); 
double f; 
double fnew; 
vector<double> g; 
vector<double> gnew; 
vector<double> xnew; 
int main() 
{ 
//Setup Timer 
clock_t start; 
clock_t end; 
double duration; 
start=clock(); 
//int dude=RSM(); 
int dude=RSM2D(); 
//End Timer 
end=clock(); 
duration=((double) (end - start)) / CLOCKS_PER_SEC; 
cout<<duration<<endl; 
system("PAUSE"); 
return(0); 
} 
double Combust(vector<double> vars) 
{ 
//New Parameters to allow for correlations 
//Define Material Properties 
double pore1=0.835; 
double pore2=vars[1]; 
double dpore1=.00029; 
double dpore2=vars[0]; 
double Omega1=0.8; 
double Omega2=0.8; 
double srho=510; 
double sCp=824; 
//Export Proerties to fill 
ofstream fid("Properties2.txt"); 
fid<<pore1<<endl; 
fid<<pore2<<endl; 
fid<<dpore1<<endl; 
fid<<dpore2<<endl; 
fid<<Omega1<<endl; 
fid<<Omega2<<endl; 
fid<<srho<<endl; 
fid<<sCp<<endl; 
fid.close(); 
//Define Fluid Properties 
double P=OneAtm; 
double Tburner=300; 
double u0=0.45; 
double phi=.65; 
IdealGasMix gas("drm19.cti","drm19"); 
int nsp=gas.nSpecies(); 
vector_fp x; 
x.resize(nsp); 
for(int k=0;k<nsp;k++) 


 { 
if(k==gas.speciesIndex("CH4")) 
{ 
x[k]=1.0; 
} 
else if(k==gas.speciesIndex("O2")) 
{ 
x[k]=0.21/phi/.105; 
} 
else if(k==gas.speciesIndex("N2")) 
{ 
x[k]=0.78/phi/.105; 
} 
else if(k==gas.speciesIndex("AR")) 
{ 
x[k]=0.01/phi/.105; 
} 
else 
{ 
x[k]=0.0; 
} 
} 
gas.setState_TPX(Tburner,P,DATA_PTR(x)); 
double rhoin=gas.density(); 
double *yin=new double[nsp]; 
gas.getMassFractions(yin); 
equilibrate(gas,"HP"); 
double rhoout=gas.density(); 
double Tad=gas.temperature(); 
double *yeq=new double[nsp]; 
gas.getMassFractions(yeq); 
//Create the grid 
double *grid=new double[301]; 
double dz1=0.0006; 
double dz2=0.00005; 
double dz3=0.00041; 
for (int i=0;i<301;i++) 
{ 
if (i<=50) 
{ 
grid[i]=((double)i)*dz1; 
} 
else if (i<=250) 
{ 
grid[i]=.03+((double)(i-50))*dz2; 
} 
else 
{ 
grid[i]=.04+((double)(i-250))*dz3; 
} 
} 
//Create the flow object 
AxiStagnFlow flow(&gas); 
flow.setupGrid(301,grid); 
Transport* tr=newTransportMgr("Mix",&gas); 
flow.setTransport(*tr); 
flow.setKinetics(gas); 
flow.setPressure(P); 
//Create the inlet 
Inlet1D inlet; 
inlet.setMoleFractions(DATA_PTR(x)); 
double mdot=u0*rhoin; 
inlet.setMdot(mdot); 
inlet.setTemperature(Tburner); 
//Create the outlet 
Outlet1D outlet; 



 //Create the flame object 
vector<Domain1D*> domains; 
domains.push_back(&inlet); 
domains.push_back(&flow); 
domains.push_back(&outlet); 
Sim1D flame(domains); 
//Build initial guess 
vector_fp locs; 
vector_fp value; 
double z1=0.55; 
double z2=0.62; 
double uout=inlet.mdot()/rhoout; 
//Velocity Profile 
locs.resize(2); 
value.resize(2); 
locs[0]=0; 
locs[1]=1; 
value[0]=u0; 
value[1]=uout; 
flame.setInitialGuess("u",locs,value); 
//Species Profiles 
locs.resize(3); 
value.resize(3); 
locs[0]=0; 
locs[1]=z1; 
locs[2]=1; 
for (int i=0;i<nsp;i++) 
{ 
value[0]=yin[i]; 
value[1]=yeq[i]; 
value[2]=yeq[i]; 
flame.setInitialGuess(gas.speciesName(i),locs,value); 
} 
//Temperature Profile 
locs.resize(4); 
value.resize(4); 
locs[0]=0; 
locs[1]=z1; 
locs[2]=z2; 
locs[3]=1; 
value[0]=Tburner; 
value[1]=Tburner; 
value[2]=2000; 
value[3]=Tad; 
flame.setInitialGuess("T",locs,value); 
//Reset Inlet 
inlet.setMoleFractions(DATA_PTR(x)); 
inlet.setMdot(mdot); 
inlet.setTemperature(Tburner); 
//Set solver parameters 
int loglevel=1; 
bool refine_grid=false; 
double rtolSS=1.0e-4; 
double atolSS=1.0e-9; 
double rtolTS=1.0e-4; 
double atolTS=1.0e-9; 
flow.setTolerancesSS(rtolSS,atolSS); 
flow.setTolerancesTS(rtolTS,atolTS); 
double SSJacAge=5; 
double TSJacAge=10; 
flame.setJacAge(SSJacAge,TSJacAge); 
//Solve and Save 
flame.solve(loglevel,refine_grid); 
refine_grid=true; 
int flowdomain=1; 
double ratio=4; 
double slope=0.4; 
double curve=0.4; 


 double prune=0.001; 
flame.setRefineCriteria(flowdomain,ratio,slope,curve,prune); 
flow.solveEnergyEqn(); 
flame.solve(loglevel,refine_grid); 
flame.save("gradienttest.xml","run","solution with energy equation"); 
flame.writeStats(); 
ifstream in("Tout.txt"); 
double input; 
in>>input; 
double Tout=input; 
double eff=pow(Tout,4)/pow(Tad,4); 
return(eff); 
} 
int RSM() 
{ 
vector<double> x0(1); 
vector<double> xprev(1); 
x0[0]=1.52; 
xprev[0]=x0[0]; 
//Selecting starting points 
vector<double> points(5); 
points[0]=1.37; 
points[1]=1.4075; 
points[2]=1.445; 
points[3]=1.4825; 
points[4]=1.52; 
double alpha=0.0375; 
//Initialize quantities 
double length=points.size(); 
vector<double> f(length); 
double diff1=10000000000000000000; 
double diff2=10000000000000000000; 
int minloc=0; 
int maxloc=4; 
double lcon=-100000000000000; //Left Constraint 
double rcon=1000000000000000; //Right Constraint 
int count=0; //Added to control number of shrinks 
double error; 
f[0]=fun2(points[0]); 
cout<<endl<<endl<<points[0]<<" "<<f[0]<<endl<<endl; 
//Loop until the model and function have same value or the model points solution is the 
same 
while ((diff1>0.0001) & (diff2>0.0001)) 
{ 
for (int i=1;i<length;i++) 
{ 
f[i]=fun2(points[i]); 
cout<<endl<<endl<<points[i]<<" "<<f[i]<<endl<<endl; 
} 
//Perform Least Squares fit 
double n=length; 
int k=1; //This is equal to the number of variables; 
double p=2*k+1; //Number of regressor variable. Need to change for higher order. 
double **X; 
X=new double* [n]; 
for (int i=0;i<n;i++) 
{ 
*(X+i)=new double[p]; 
} 
for (int i=0;i<n;i++) 
{ 
for (int c=0;c<k;c++) 
{ 
X[i][c]=1; 
X[i][c+1]=points[i]; 
X[i][c+k+1]=pow(points[i],2); 
} 


 } 
double **Xt; 
Xt=new double* [p]; 
for (int i=0;i<p;i++) 
{ 
*(Xt+i)=new double[n]; 
} 
Xt=MT(X,n,p); 
double **A; 
A=new double* [p]; 
for (int i=0;i<p;i++) 
{ 
*(A+i)=new double[p]; 
} 
A=MM(Xt,X,p,n,p); 
vector<double> c(p); 
c=Mv(Xt,f,p,n); 
vector<double> b(p); 
b=LUSolve(A,c); 
//Perform Newtons Method 
double g; 
double H; 
double fmin=newfun(x0,b,g,H); 
if (H<0) 
{ 
double g1; 
double H1; 
vector<double> p(1); 
p[0]=points[minloc]; 
double b1=newfun(p,b,g1,H1); 
double g2; 
double H2; 
p[0]=points[maxloc]; 
double b2=newfun(p,b,g2,H2); 
if (b1<=b2) 
{ 
fmin=b1; 
x0[0]=points[minloc]; 
} 
else 
{ 
fmin=b2; 
x0[0]=points[maxloc]; 
} 
} 
else 
{ 
double d=-g/H; 
double temp; //Added to stop solver from going beyond model range. 
temp=x0[0]+d; 
if (temp<points[minloc]) 
{ 
x0[0]=points[minloc]; 
} 
else if (temp>points[maxloc]) 
{ 
x0[0]=points[maxloc]; 
} 
else 
{ 
x0[0]=x0[0]+d; 
if (count<3) 
{ 
alpha=alpha/2;//Shrinks if convex and inside box 
count=count+1; 
} 
else 
{ 
xprev[0]=x0[0]; 


 } 
} 
fmin=newfun(x0,b,g,H); 
} 
cout<<x0[0]<<endl<<endl; 
double fnew=fun2(x0[0]); 
int ignore=0; 
if (fnew>f[0]) 
{ 
x0[0]=points[0]; 
fnew=f[0]; 
xprev[0]=x0[0]; 
ignore=1; 
} 
//Error Stuff 
if (ignore==0) 
{ 
vector<double> xm(p); 
xm[0]=1; 
xm[1]=x0[0]; 
xm[2]=pow(x0[0],2); 
vector<double> yhat(n); 
yhat=Mv(X,b,n,p); 
vector<double> fsurf(n); 
for (int i=0;i<n;i++) 
{ 
vector<double> node(1); 
node[0]=points[0]; 
fsurf[i]=pow((f[i]-yhat[i]),2)/(n-p); 
} 
double temp=sum(fsurf); 
double s=sqrt(temp); 
vector<double> error1(p); 
error1=LUSolve(A,xm); 
error=inprod(xm,error1); 
error=s*2.919986*sqrt(error); //The number is for 90% confidence from 
students t 
} 
cout<<endl<<endl<<x0[0]<<" "<<fnew<<endl<<endl; 
//Update loop ending parameters 
diff1=sqrt(pow((fnew-fmin)/fnew,2)); 
vector<double> xtemp(1); 
xtemp[0]=x0[0]-xprev[0]; 
xprev[0]=x0[0]; 
diff2=norm(xtemp); 
cout<<diff1<<" "<<diff2<<endl; 
f[0]=fnew; 
//Update points 
points[0]=x0[0]; 
if (x0[0]<(lcon+2*alpha)) 
{ 
if (x0[0]<=(lcon+alpha)) 
{ 
if (x0[0]<=lcon) 
{ 
if (count<3) 
{ 
alpha=alpha/2;//Shrink if on min edge 
count=count+1; 
} 
points[0]=lcon; 
points[1]=points[0]+alpha; 
points[2]=points[1]+alpha; 
points[3]=points[2]+alpha; 
points[4]=points[3]+alpha; 
minloc=0; 


 maxloc=4; 
} 
else 
{ 
points[1]=lcon; 
points[2]=points[0]+alpha; 
points[3]=points[2]+alpha; 
points[4]=points[3]+alpha; 
minloc=1; 
maxloc=4; 
} 
} 
else 
{ 
points[1]=lcon; 
points[2]=points[0]-alpha; 
points[3]=points[0]+alpha; 
points[4]=points[3]+alpha; 
minloc=1; 
maxloc=4; 
} 
} 
else if (x0[0]>(rcon-2*alpha)) 
{ 
if (x0[0]>=(rcon-alpha)) 
{ 
if (x0[0]>=rcon) 
{ 
if (count<3) 
{ 
alpha=alpha/2;//Shrink if on max edge 
count=count+1; 
} 
points[0]=rcon; 
points[1]=points[0]-alpha; 
points[2]=points[1]-alpha; 
points[3]=points[2]-alpha; 
points[4]=points[3]-alpha; 
minloc=4; 
maxloc=0; 
} 
else 
{ 
points[1]=rcon; 
points[2]=points[0]-alpha; 
points[3]=points[2]-alpha; 
points[4]=points[3]-alpha; 
minloc=4; 
maxloc=1; 
} 
} 
else 
{ 
points[1]=rcon; 
points[2]=points[0]+alpha; 
points[3]=points[0]-alpha; 
points[4]=points[3]-alpha; 
minloc=4; 
maxloc=1; 
} 
} 
else 
{ 
points[1]=points[0]-alpha; 
points[2]=points[1]-alpha; 
points[3]=points[0]+alpha; 
points[4]=points[3]+alpha; 
minloc=2; 
maxloc=4; 
} 
} 


 cout<<"x*="<<x0[0]<<" with fmax="<<f[0]<<"+/-"<<error<<endl; 
return(0); 
} 
int RSM2D() 
{ 
vector<double> x0(2); 
vector<double> xprev(2); 
x0[0]=1.445; 
x0[1]=0.88; 
xprev[0]=x0[0]; 
xprev[1]-x0[1]; 
//Selecting starting points 
double **points; 
points=new double* [9]; 
for (int i=0;i<9;i++) 
{ 
*(points+i)=new double[2]; 
} 
points[0][0]=1.37; 
points[1][0]=1.37; 
points[2][0]=1.37; 
points[3][0]=1.445; 
points[4][0]=1.445; 
points[5][0]=1.445; 
points[6][0]=1.52; 
points[7][0]=1.52; 
points[8][0]=1.52; 
points[0][1]=.87; 
points[1][1]=.88; 
points[2][1]=.89; 
points[3][1]=.87; 
points[4][1]=.88; 
points[5][1]=.89; 
points[6][1]=.87; 
points[7][1]=.88; 
points[8][1]=.89; 
double alpha1=0.075; 
double alpha2=0.01; 
//Initialize quantities 
double length=9; 
vector<double> f(length); 
double diff1=10000000000000000000; 
double diff2=10000000000000000000; 
double lcon=0.69; 
double rcon=1.52; 
double bcon=0.865; 
double tcon=0.95; 
double error; 
int flagl=0; 
int flagr=0; 
int flagb=0; 
int flagt=0; 
int count1=0; 
int count2=0; 
int nshrinks=4; 
int noshrinkh=0; 
int noshrinkv=0; 
int start=0; 
f[0]=fun2D(points[0][0],points[0][1]); 
//Loop until the model points solution is the same 
while (diff2>0.0001) 
{ 
for (int i=1;i<length;i++) 
{ 
f[i]=fun2D(points[i][0],points[i][1]); 
} 
//Perform Least Squares fit 


 double n=length; 
int k=2; 
double p=2*k+2; 
double **X; 
X=new double* [n]; 
for (int i=0;i<n;i++) 
{ 
*(X+i)=new double[p]; 
} 
for (int i=0;i<n;i++) 
{ 
X[i][0]=1; 
for (int c=0;c<k;c++) 
{ 
X[i][c+1]=points[i][c]; 
X[i][c+k+1]=pow(points[i][c],2); 
} 
for (int c=0;c<k-1;c++) 
{ 
for (int j=c;j<k;j++) 
{ 
X[i][c+2*k+1]=points[i][c]*points[i][j]; 
} 
} 
} 
double **Xt; 
Xt=new double* [p]; 
for (int i=0;i<p;i++) 
{ 
*(Xt+i)=new double[n]; 
} 
Xt=MT(X,n,p); 
double **A; 
A=new double* [p]; 
for (int i=0;i<p;i++) 
{ 
*(A+i)=new double[p]; 
} 
A=MM(Xt,X,p,n,p); 
vector<double> c(p); 
c=Mv(Xt,f,p,n); 
vector<double> b(p); 
b=LUSolve(A,c); 
//Perform Newtons Method 
vector<double> grad(2); 
double **Hess; 
Hess=new double* [2]; 
for (int i=0;i<2;i++) 
{ 
*(Hess+i)=new double[2]; 
} 
double fmin=newfun2D(x0,b,grad,Hess); 
vector<double> d(2); 
vector<double> gneg(2); 
gneg[0]=-grad[0]; 
gneg[1]=-grad[1]; 
double **Hes; 
Hes=new double* [2]; 
for (int i=0;i<2;i++) 
{ 
*(Hes+i)=new double[2]; 
} 
Hes[0][0]=Hess[0][0]; 
Hes[0][1]=Hess[0][1]; 
Hes[1][0]=Hess[1][0]; 
Hes[1][1]=Hess[1][1]; 


 d=LUSolve(Hes,gneg); 
//New code added to deal with indefinite Hessians. Gaurentees descent. 
double LS=inprod(gneg,d); 
if (LS<0) 
{ 
vector<double> dprime(2); 
dprime[0]=sqrt(pow(gneg[0],2)); 
dprime[1]=sqrt(pow(gneg[1],2)); 
double alpha1prime=min(alpha1,min(rcon-x0[0],x0[0]-lcon)); 
double alpha2prime=min(alpha2,min(tcon-x0[1],x0[1]-bcon)); 
if ((alpha1prime==rcon-x0[0])&(sign(gneg[0])==-1)) 
{ 
alpha1prime=alpha1; 
} 
else if ((alpha1prime==x0[0]-lcon)&(sign(gneg[0])==1)) 
{ 
alpha1prime=alpha1; 
} 
if ((alpha2prime==tcon-x0[1])&(sign(gneg[1])==-1)) 
{ 
alpha2prime=alpha2; 
} 
else if ((alpha2prime==x0[1]-bcon)&(sign(gneg[1])==1)) 
{ 
alpha2prime=alpha2; 
} 
vector<double> dnew(2); 
double gamma1=atan(alpha2prime/alpha1prime); 
double gamma2=atan(dprime[1]/dprime[0]); 
if (gamma2>=gamma1) 
{ 
dnew[0]=dprime[0]*alpha2prime/dprime[1]; 
dnew[1]=alpha2prime; 
d[0]=sign(gneg[0])*dnew[0]; 
d[1]=sign(gneg[1])*dnew[1]; 
} 
else 
{ 
dnew[0]=alpha1prime; 
dnew[1]=dprime[1]*alpha1prime/dprime[0]; 
d[0]=sign(gneg[0])*dnew[0]; 
d[1]=sign(gneg[1])*dnew[1]; 
} 
} 
else 
{ 
vector<double> dprime(2); 
dprime[0]=sqrt(pow(d[0],2)); 
dprime[1]=sqrt(pow(d[1],2)); 
double alpha1prime=min(alpha1,min(rcon-x0[0],x0[0]-lcon)); 
double alpha2prime=min(alpha2,min(tcon-x0[1],x0[1]-bcon)); 
if ((alpha1prime==rcon-x0[0])&(sign(d[0])==-1)) 
{ 
alpha1prime=alpha1; 
} 
else if ((alpha1prime==x0[0]-lcon)&(sign(d[0])==1)) 
{ 
alpha1prime=alpha1; 
} 
if ((alpha2prime==tcon-x0[1])&(sign(d[1])==-1)) 
{ 
alpha2prime=alpha2; 
} 
else if ((alpha2prime==x0[1]-bcon)&(sign(d[1])==1)) 
{ 
alpha2prime=alpha2; 
} 


 vector<double> dnew(2); 
if ((dprime[0]>alpha1prime)|(dprime[1]>alpha2prime)) 
{ 
double gamma1=atan(alpha2prime/alpha1prime); 
double gamma2=atan(dprime[1]/dprime[0]); 
if (gamma2>=gamma1) 
{ 
dnew[0]=dprime[0]*alpha2prime/dprime[1]; 
dnew[1]=alpha2prime; 
d[0]=sign(d[0])*dnew[0]; 
d[1]=sign(d[1])*dnew[1]; 
} 
else 
{ 
dnew[0]=alpha1prime; 
dnew[1]=dprime[1]*alpha1prime/dprime[0]; 
d[0]=sign(d[0])*dnew[0]; 
d[1]=sign(d[1])*dnew[1]; 
} 
} 
else 
{ 
if (count1<nshrinks) 
{ 
if (count2<nshrinks) 
{ 
alpha1=alpha1/2; 
alpha2=alpha2/2; 
count1=count1+1; 
count2=count2+1; 
} 
else 
{ 
alpha1=alpha1/2; 
count1=count1+1; 
} 
} 
else if (count2<nshrinks) 
{ 
alpha2/alpha2/2; 
count2=count2+1; 
} 
else 
{ 
xprev[0]=x0[0]; 
xprev[1]=x0[1]; 
} 
} 
} 
//Update point 
x0[0]=x0[0]+d[0]; 
x0[1]=x0[1]+d[1]; 
fmin=newfun2D(x0,b,grad,Hess); 
double fnew=fun2D(x0[0],x0[1]); 
cout<<endl<<endl<<x0[0]<<" "<<x0[1]<<" "<<fnew<<endl<<endl; 
//If new point is worse than previous stop 
int ignore=0; 
if (fnew>f[0]) 
{ 
if (start==1) 
{ 
x0[0]=points[0][0]; 
x0[1]=points[0][1]; 
fnew=f[0]; 
flagt=0; 
flagb=0; 
flagl=0; 
flagr=0; 
xprev[0]=x0[0]; 


 xprev[1]=x0[1]; 
diff2=0; 
ignore=1; 
} 
else 
{ 
start=1; 
} 
} 
//Error Calc 
if (ignore==0) 
{ 
vector<double> xm(p); 
xm[0]=1; 
xm[1]=x0[0]; 
xm[2]=x0[1]; 
xm[3]=pow(x0[0],2); 
xm[4]=pow(x0[1],2); 
xm[5]=x0[0]*x0[1]; 
vector<double> yhat(n); 
yhat=Mv(X,b,n,p); 
vector<double> fsurf(n); 
for (int i=0;i<n;i++) 
{ 
fsurf[i]=pow((f[i]-yhat[i]),2)/(n-p); 
} 
double temp=sum(fsurf); 
double s=sqrt(temp); 
vector<double> error1(p); 
error1=LUSolve(A,xm); 
error=inprod(xm,error1); 
error=s*2.919986*sqrt(error); //The number is for 90% confidence from 
students t 
} 
//Update loop ending parameters 
diff1=sqrt(pow((fnew-fmin)/fnew,2)); 
vector<double> xtemp(2); 
xtemp[0]=x0[0]-xprev[0]; 
xtemp[1]=x0[1]-xprev[1]; 
xprev[0]=x0[0]; 
xprev[1]=x0[1]; 
diff2=norm(xtemp); 
cout<<diff1<<" "<<diff2<<endl; 
//Check if GRGM needs to be used 
if (diff2<0.0001) 
{ 
if ((flagl==1)|(flagr==1)) 
{ 
fnew=GRGM(x0,fnew,alpha2,bcon,tcon,1); 
if (count2<nshrinks) 
{ 
count2=count2+1; 
alpha2=alpha2/2; 
if (noshrinkh==0) 
{ 
alpha1=alpha1*2; 
noshrinkh=1; 
} 
} 
else 
{ 
xprev[0]=x0[0]; 
xprev[1]=x0[1]; 
} 
} 
else if ((flagb==1)|(flagt==1)) 
{ 
fnew=GRGM(x0,fnew,alpha1,lcon,rcon,0); 


 if (count1<nshrinks) 
{ 
count1=count1+1; 
alpha1=alpha1/2; 
if (noshrinkv==0) 
{ 
alpha2=alpha2*2; 
noshrinkv=1; 
} 
} 
else 
{ 
xprev[0]=x0[0]; 
xprev[1]=x0[1]; 
} 
} 
xtemp[0]=x0[0]-xprev[0]; 
xtemp[1]=x0[1]-xprev[1]; 
xprev[0]=x0[0]; 
xprev[1]=x0[1]; 
diff2=norm(xtemp); 
cout<<diff1<<" "<<diff2<<endl; 
//system("PAUSE"); 
} 
f[0]=fnew; 
//Update points 
points[0][0]=x0[0]; 
points[0][1]=x0[1]; 
if (x0[0]<(lcon+alpha1)) 
{ 
if (x0[0]<=lcon) 
{ 
if (count1<nshrinks) 
{ 
alpha1=alpha1/2;//Shrink if on left edge. 
count1=count1+1; 
flagl=1; 
} 
else 
{ 
flagl=0; 
} 
points[0][0]=lcon; 
points[1][0]=lcon; 
points[2][0]=lcon; 
points[3][0]=points[0][0]+alpha1; 
points[4][0]=points[0][0]+alpha1; 
points[5][0]=points[0][0]+alpha1; 
points[6][0]=points[3][0]+alpha1; 
points[7][0]=points[3][0]+alpha1; 
points[8][0]=points[3][0]+alpha1; 
flagr=0; 
flagb=0; 
flagt=0; 
} 
else 
{ 
points[1][0]=points[0][0]; 
points[2][0]=points[0][0]; 
points[3][0]=lcon; 
points[4][0]=lcon; 
points[5][0]=lcon; 
points[6][0]=points[0][0]+alpha1; 
points[7][0]=points[0][0]+alpha1; 
points[8][0]=points[0][0]+alpha1; 
flagl=0; 
flagr=0; 
flagb=0; 
flagt=0; 


 noshrinkh=0; 
} 
} 
else if (x0[0]>(rcon-alpha1)) 
{ 
if (x0[0]>=rcon) 
{ 
if (count1<nshrinks) 
{ 
alpha1=alpha1/2;//Shrink if on right edge. 
count1=count1+1; 
flagr=1; 
} 
else 
{ 
flagr=0; 
} 
points[0][0]=rcon; 
points[1][0]=rcon; 
points[2][0]=rcon; 
points[3][0]=points[0][0]-alpha1; 
points[4][0]=points[0][0]-alpha1; 
points[5][0]=points[0][0]-alpha1; 
points[6][0]=points[3][0]-alpha1; 
points[7][0]=points[3][0]-alpha1; 
points[8][0]=points[3][0]-alpha1; 
flagl=0; 
flagb=0; 
flagt=0; 
} 
else 
{ 
points[1][0]=points[0][0]; 
points[2][0]=points[0][0]; 
points[3][0]=rcon; 
points[4][0]=rcon; 
points[5][0]=rcon; 
points[6][0]=points[0][0]-alpha1; 
points[7][0]=points[0][0]-alpha1; 
points[8][0]=points[0][0]-alpha1; 
flagl=0; 
flagr=0; 
flagb=0; 
flagt=0; 
noshrinkh=0; 
} 
} 
else 
{ 
points[1][0]=points[0][0]; 
points[2][0]=points[0][0]; 
points[3][0]=points[0][0]+alpha1; 
points[4][0]=points[0][0]+alpha1; 
points[5][0]=points[0][0]+alpha1; 
points[6][0]=points[0][0]-alpha1; 
points[7][0]=points[0][0]-alpha1; 
points[8][0]=points[0][0]-alpha1; 
flagl=0; 
flagr=0; 
flagb=0; 
flagt=0; 
noshrinkh=0; 
} 
if (x0[1]<(bcon+alpha2)) 
{ 
if (x0[1]<=bcon) 
{ 
if (count2<nshrinks) 
{ 
alpha2=alpha2/2;//Shrink if on bottom edge. 


 count2=count2+1; 
flagb=1; 
} 
else 
{ 
flagb=0; 
} 
points[0][1]=bcon; 
points[1][1]=points[0][1]+alpha2; 
points[2][1]=points[1][1]+alpha2; 
points[3][1]=bcon; 
points[4][1]=points[0][1]+alpha2; 
points[5][1]=points[1][1]+alpha2; 
points[6][1]=bcon; 
points[7][1]=points[0][1]+alpha2; 
points[8][1]=points[1][1]+alpha2; 
flagl=0; 
flagr=0; 
flagt=0; 
} 
else 
{ 
points[1][1]=points[0][1]+alpha2; 
points[2][1]=bcon; 
points[3][1]=points[0][1]; 
points[4][1]=points[0][1]+alpha2; 
points[5][1]=bcon; 
points[6][1]=points[0][1]; 
points[7][1]=points[0][1]+alpha2; 
points[8][1]=bcon; 
flagl=0; 
flagr=0; 
flagb=0; 
flagt=0; 
noshrinkv=0; 
} 
} 
else if (x0[1]>(tcon-alpha2)) 
{ 
if (x0[1]>=tcon) 
{ 
if (count2<nshrinks) 
{ 
alpha2=alpha2/2;//Shrink if on top edge. 
count2=count2+1; 
flagt=1; 
} 
else 
{ 
flagt=0; 
} 
points[0][1]=tcon; 
points[1][1]=points[0][1]-alpha2; 
points[2][1]=points[1][1]-alpha2; 
points[3][1]=tcon; 
points[4][1]=points[0][1]-alpha2; 
points[5][1]=points[1][1]-alpha2; 
points[6][1]=tcon; 
points[7][1]=points[0][1]-alpha2; 
points[8][1]=points[1][1]-alpha2; 
flagl=0; 
flagr=0; 
flagb=0; 
} 
else 
{ 
points[1][1]=points[0][1]-alpha2; 
points[2][1]=tcon; 
points[3][1]=points[0][1]; 
points[4][1]=points[0][1]-alpha2; 
points[5][1]=tcon; 


 points[6][1]=points[0][1]; 
points[7][1]=points[0][1]-alpha2; 
points[8][1]=tcon; 
flagl=0; 
flagr=0; 
flagb=0; 
flagt=0; 
noshrinkv=0; 
} 
} 
else 
{ 
points[1][1]=points[0][1]+alpha2; 
points[2][1]=points[0][1]-alpha2; 
points[3][1]=points[0][1]; 
points[4][1]=points[0][1]+alpha2; 
points[5][1]=points[0][1]-alpha2; 
points[6][1]=points[0][1]; 
points[7][1]=points[0][1]+alpha2; 
points[8][1]=points[0][1]-alpha2; 
flagl=0; 
flagr=0; 
flagb=0; 
flagt=0; 
noshrinkv=0; 
} 
} 
cout<<"x*=("<<x0[0]<<","<<x0[1]<<") with fmax="<<f[0]<<"+/-"<<error<<endl; 
return(0); 
} 
double GRGM(vector<double> &x0, double f0, double alpha, double lcon, double rcon, int dim) 
{ 
//Selecting starting points 
alpha=alpha/2; 
int stop=0; 
double minloc; 
double maxloc; 
int adim; 
if (dim==1) 
{ 
adim=0; 
} 
else 
{ 
adim=1; 
} 
vector<double> points(5); 
points[0]=x0[dim]; 
//checks which piece of the x vector we are altering and then chooses points 
if (x0[dim]<(lcon+2*alpha)) 
{ 
if (x0[dim]<=(lcon+alpha)) 
{ 
if (x0[dim]<=lcon) 
{ 
stop=1; 
} 
else 
{ 
points[1]=lcon; 
points[2]=points[0]+alpha; 
points[3]=points[2]+alpha; 
points[4]=points[3]+alpha; 
minloc=1; 
maxloc=4; 
} 
} 
else 
{ 


 points[1]=lcon; 
points[2]=points[0]-alpha; 
points[3]=points[0]+alpha; 
points[4]=points[3]+alpha; 
minloc=1; 
maxloc=4; 
} 
} 
else if (x0[dim]>(rcon-2*alpha)) 
{ 
if (x0[dim]>=(rcon-alpha)) 
{ 
if (x0[dim]>=rcon) 
{ 
stop=1; 
} 
else 
{ 
points[1]=rcon; 
points[2]=points[0]-alpha; 
points[3]=points[2]-alpha; 
points[4]=points[3]-alpha; 
minloc=4; 
maxloc=1; 
} 
} 
else 
{ 
points[1]=rcon; 
points[2]=points[0]+alpha; 
points[3]=points[0]-alpha; 
points[4]=points[3]-alpha; 
minloc=4; 
maxloc=1; 
} 
} 
else 
{ 
points[1]=points[0]-alpha; 
points[2]=points[1]-alpha; 
points[3]=points[0]+alpha; 
points[4]=points[3]+alpha; 
minloc=2; 
maxloc=4; 
} 
//Initialize quantities 
double length=points.size(); 
double fnew; 
vector<double> f(length); 
f[0]=f0; 
//Loop until the model and function have same value or the model points solution is the 
same 
while (stop==0) 
{ 
for (int i=1;i<length;i++) 
{ 
if (adim==0) 
{ 
f[i]=fun2D(x0[adim],points[i]); 
} 
else 
{ 
f[i]=fun2D(points[i],x0[adim]); 
} 
} 
//Perform Least Squares fit 
double n=length; 
int k=1; //This is equal to the number of variables; 
double p=2*k+1; //Number of regressor variable. Need to change for higher order. 



 double **X; 
X=new double* [n]; 
for (int i=0;i<n;i++) 
{ 
*(X+i)=new double[p]; 
} 
for (int i=0;i<n;i++) 
{ 
//Need to add here for higher dimensions 
for (int c=0;c<k;c++) 
{ 
X[i][c]=1; 
X[i][c+1]=points[i]; 
X[i][c+k+1]=pow(points[i],2); 
} 
} 
double **Xt; 
Xt=new double* [p]; 
for (int i=0;i<p;i++) 
{ 
*(Xt+i)=new double[n]; 
} 
Xt=MT(X,n,p); 
double **A; 
A=new double* [p]; 
for (int i=0;i<p;i++) 
{ 
*(A+i)=new double[p]; 
} 
A=MM(Xt,X,p,n,p); 
vector<double> c(p); 
c=Mv(Xt,f,p,n); 
vector<double> b(p); 
b=LUSolve(A,c); 
//Perform Newtons Method 
double g; 
double H; 
double fmin=newfun(x0,b,g,H); 
if (H<0) 
{ 
double g1; 
double H1; 
vector<double> p(1); 
p[0]=points[minloc]; 
double b1=newfun(p,b,g1,H1); 
double g2; 
double H2; 
p[0]=points[maxloc]; 
double b2=newfun(p,b,g2,H2); 
if (b1<=b2) 
{ 
fmin=b1; 
x0[dim]=points[minloc]; 
} 
else 
{ 
fmin=b2; 
x0[dim]=points[maxloc]; 
} 
} 
else 
{ 
double d=-g/H; 
double temp; //Added to stop solver from going beyond model range. 
temp=x0[dim]+d; 
if (temp<points[minloc]) 
{ 
x0[dim]=points[minloc]; 


 } 
else if (temp>points[maxloc]) 
{ 
x0[dim]=points[maxloc]; 
} 
else 
{ 
x0[dim]=x0[dim]+d; 
stop=1; 
} 
} 
if (adim==0) 
{ 
fnew=fun2D(x0[adim],x0[dim]); 
} 
else 
{ 
fnew=fun2D(x0[dim],x0[adim]); 
} 
cout<<endl<<endl<<x0[dim]<<" "<<fnew<<endl<<endl; 
//system("PAUSE"); 
if (fnew>f[0]) 
{ 
x0[dim]=points[0]; 
fnew=f[0]; 
stop=1; 
} 
f[0]=fnew; 
//Update points 
points[0]=x0[dim]; 
//Added consraints back in and generalized April 14, 2010 
if (x0[dim]<(lcon+2*alpha)) 
{ 
if (x0[dim]<=(lcon+alpha)) 
{ 
if (x0[dim]<=lcon) 
{ 
stop=1; 
} 
else 
{ 
points[1]=lcon; 
points[2]=points[0]+alpha; 
points[3]=points[2]+alpha; 
points[4]=points[3]+alpha; 
minloc=1; 
maxloc=4; 
} 
} 
else 
{ 
points[1]=lcon; 
points[2]=points[0]-alpha; 
points[3]=points[0]+alpha; 
points[4]=points[3]+alpha; 
minloc=1; 
maxloc=4; 
} 
} 
else if (x0[dim]>(rcon-2*alpha)) 
{ 
if (x0[dim]>=(rcon-alpha)) 
{ 
if (x0[dim]>=rcon) 
{ 
stop=1; 
} 
else 
{ 
points[1]=rcon; 


 points[2]=points[0]-alpha; 
points[3]=points[2]-alpha; 
points[4]=points[3]-alpha; 
minloc=4; 
maxloc=1; 
} 
} 
else 
{ 
points[1]=rcon; 
points[2]=points[0]+alpha; 
points[3]=points[0]-alpha; 
points[4]=points[3]-alpha; 
minloc=4; 
maxloc=1; 
} 
} 
else 
{ 
points[1]=points[0]-alpha; 
points[2]=points[1]-alpha; 
points[3]=points[0]+alpha; 
points[4]=points[3]+alpha; 
minloc=2; 
maxloc=4; 
} 
} 
return(fnew); 
} 
double norm(vector<double> x) 
//finds norms of vectors 
{ 
double length=x.size(); 
double sum=0; 
for (int i=0;i<length;i++) 
{ 
sum+=x[i]*x[i]; 
} 
double norm=sqrt(sum); 
return(norm); 
} 
double inprod(vector<double>x, vector<double>y) 
//inner product of vectors 
{ 
double length=x.size(); 
double ans=0; 
for(int i=0;i<length;i++) 
{ 
ans+=x[i]*y[i]; 
} 
return(ans); 
} 
vector<double> Mv(double **M,vector<double> v,double row, double col) 
//matrix times a vector 
{ 
vector<double> ans(row,0); 
for(int i=0;i<row;i++) 
{ 
for(int j=0;j<col;j++) 
{ 
ans[i]=ans[i]+M[i][j]*v[j]; 
} 
} 
return(ans); 
} 
double **MM(double **M1, double **M2,double a, double b, double c) 
//matrix times matrix 
{ 


 double **ans; 
ans=new double* [a]; 
for(int i=0;i<a;i++) 
{ 
*(ans+i)=new double[c]; 
} 
for(int i=0;i<a;i++) 
{ 
for(int j=0;j<c;j++) 
{ 
ans[i][j]=0; 
} 
} 
for(int k=0;k<a;k++) 
{ 
for(int i=0;i<c;i++) 
{ 
for(int j=0;j<b;j++) 
{ 
ans[k][i]=ans[k][i]+M1[k][j]*M2[j][i]; 
} 
} 
} 
return(ans); 
} 
vector<double> LUSolve(double **M,vector<double> v) 
//LU decomposion solver 
{ 
int length=v.size(); 
for(int i=0;i<length-1;i++) 
{ 
for(int j=i+1;j<length;j++) 
{ 
double m=M[j][i]/M[i][i]; 
M[j][i]=0; 
for(int k=i+1;k<length;k++) 
{ 
M[j][k]=M[j][k]-m*M[i][k]; 
} 
M[j][i]=m; 
} 
} 
for(int i=1;i<length;i++) 
{ 
for (int j=0;j<i;j++) 
{ 
v[i]=v[i]-M[i][j]*v[j]; 
} 
} 
v[length-1]=v[length-1]/M[length-1][length-1]; 
for (int i=0;i<length-1;i++) 
{ 
for (int j=length-1-i;j<length;j++) 
{ 
v[length-2-i]=v[length-2-i]-M[length-2-i][j]*v[j]; 
} 
v[length-2-i]=v[length-2-i]/M[length-2-i][length-2-i]; 
} 
return(v); 
} 
double sum(vector<double> x) 
//add vector components 
{ 
double add=0; 
double length=x.size(); 
for (int i=0;i<length;i++) 
{ 
add+=x[i]; 
} 
return (add); 


} 
double newfun(vector<double> x, vector<double> b, double &g, double &H) 
//1-D model function 
{ 
double val=b[0]+b[1]*x[0]+b[2]*pow(x[0],2); 
g=b[1]+2*b[2]*x[0]; 
H=2*b[2]; 
return(val); 
} 
double newfun2D(vector<double> x, vector<double> b, vector<double> &grad, double **&Hess) 
//2-D model function 
{ 
double val=b[0]+b[1]*x[0]+b[2]*x[1]+b[3]*pow(x[0],2)+b[4]*pow(x[1],2)+b[5]*x[0]*x[1]; 
grad[0]=b[1]+2*b[3]*x[0]+b[5]*x[1]; 
grad[1]=b[2]+2*b[4]*x[1]+b[5]*x[0]; 
Hess[0][0]=2*b[3]; 
Hess[0][1]=b[5]; 
Hess[1][0]=b[5]; 
Hess[1][1]=2*b[4]; 
return(val); 
} 
double fun2(double x) 
//combustion function evaluation for 1-D case 
{ 
vector<double> vars(1); 
vars[0]=x/1000; //Divide by 1000 for pore diameter 
double eff=-Combust(vars); //Added negative for maximization. 
return(eff); 
} 
double fun2D(double x, double y) 
//combustion function evaluation for 2-D case 
{ 
//double eff=100*pow((y-pow(x,2)),2)+pow((1-x),2); //Rosenbrock 
vector<double> vars(2); 
vars[0]=x/1000; 
vars[1]=y; 
double eff=-Combust(vars); 
return(eff); 
} 
double **MT(double **M,double n, double p) 
//matrix times its transpose 
{ 
double **ans; 
ans=new double* [p]; 
for(int i=0;i<p;i++) 
{ 
*(ans+i)=new double[n]; 
} 
for (int i=0;i<p;i++) 
{ 
for (int j=0;j<n;j++) 
{ 
ans[i][j]=M[j][i]; 
} 
} 
return(ans); 
} 
double max(double x, double y) 
//max of two numbers 
{ 
double temp1=sqrt(pow(x,2)); 
double temp2=sqrt(pow(y,2)); 
double ans; 
if (temp1>=temp2) 
{ 
ans=temp1; 
} 
else 
{ 
ans=temp2; 
} 


 return(ans); 
} 
double min(double x, double y) 
//min of two numbers 
{ 
double temp1=sqrt(pow(x,2)); 
double temp2=sqrt(pow(y,2)); 
double ans; 
if (temp1<=temp2) 
{ 
ans=temp1; 
} 
else 
{ 
ans=temp2; 
} 
return(ans); 
} 
int sign(double x) 
//sign of a number 
{ 
int ans; 
if (x>=0) 
{ 
ans=1; 
} 
else 
{ 
ans=-1; 
} 
return(ans); 

}
