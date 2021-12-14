C.4 Stflow.cpp 
/** 
* @file StFlow.cpp 
*/ 
/* 
* $Author: dggoodwin $ 
* $Revision: 1.29 $ 
* $Date: 2006/04/28 17:22:23 $ 
*/ 
// Copyright 2002 California Institute of Technology 
// turn off warnings under Windows 
#ifdef WIN32 
#pragma warning(disable:4786) 
#pragma warning(disable:4503) 
#pragma warning(disable:4267) 
#endif 
#include <stdlib.h> 
#include <time.h> 
#include <vector> 
#include <fstream> 
#include "StFlow.h" 
#include "../ArrayViewer.h" 
#include "../ctml.h" 
#include "MultiJac.h" 
#include "OneDim.cpp" 
using namespace ctml; 
using namespace std; 
namespace Cantera { 
//initialize solid temperature vector, the previous temperature profile, the previous 
mesh, 
//the heat transfer coefficient, and no adaption 
vector<double> Tw; 
vector<double> Twprev; 
vector<double> Twprev1; 
vector<double> zprev; 
vector<double> hconv; 
int adapt=0; 
//------------------- importSolution ------------------------ 
/** 
* Import a previous solution to use as an initial estimate. The 
* previous solution may have been computed using a different 
* reaction mechanism. Species in the old and new mechanisms are 
* matched by name, and any species in the new mechanism that were 
* not in the old one are set to zero. The new solution is created 
* with the same number of grid points as in the old solution. 
*/ 
void importSolution(int points, 
doublereal* oldSoln, igthermo_t& oldmech, 
int size_new, doublereal* newSoln, igthermo_t& newmech) { 
// Number of components in old and new solutions 
int nv_old = oldmech.nSpecies() + 4; 
int nv_new = newmech.nSpecies() + 4; 
if (size_new < nv_new*points) { 
throw CanteraError("importSolution", 
"new solution array must have length "+ 


 int2str(nv_new*points)); 
} 
int n, j, knew; 
string nm; 
// copy u,V,T,lambda 
for (j = 0; j < points; j++) 
for (n = 0; n < 4; n++) 
newSoln[nv_new*j + n] = oldSoln[nv_old*j + n]; 
// copy mass fractions 
int nsp0 = oldmech.nSpecies(); 
//int nsp1 = newmech.nSpecies(); 
// loop over the species in the old mechanism 
for (int k = 0; k < nsp0; k++) { 
nm = oldmech.speciesName(k); // name 
// location of this species in the new mechanism. 
// If < 0, then the species is not in the new mechanism. 
knew = newmech.speciesIndex(nm); 
// copy this species from the old to the new solution vectors 
if (knew >= 0) { 
for (j = 0; j < points; j++) { 
newSoln[nv_new*j + 4 + knew] = oldSoln[nv_old*j + 4 + k]; 
} 
} 
} 
// normalize mass fractions 
for (j = 0; j < points; j++) { 
newmech.setMassFractions(&newSoln[nv_new*j + 4]); 
newmech.getMassFractions(&newSoln[nv_new*j + 4]); 
} 
} 
static void st_drawline() { 
writelog("\n-------------------------------------" 
"------------------------------------------"); 
} 
StFlow::StFlow(igthermo_t* ph, int nsp, int points) : 
Domain1D(nsp+4, points), 
m_inlet_u(0.0), 
m_inlet_V(0.0), 
m_inlet_T(-1.0), 
m_surface_T(-1.0), 
m_press(-1.0), 
m_nsp(nsp), 
m_thermo(0), 
m_kin(0), 
m_trans(0), 
m_jac(0), 
m_ok(false), 
m_do_soret(false), 
m_transport_option(-1), 
m_efctr(0.0) 
{ 
m_type = cFlowType; 
m_points = points; 
m_thermo = ph; 
if (ph == 0) return; // used to create a dummy object 
int nsp2 = m_thermo->nSpecies(); 
if (nsp2 != m_nsp) { 


 m_nsp = nsp2; 
Domain1D::resize(m_nsp+4, points); 
} 
// make a local copy of the species molecular weight vector 
m_wt = m_thermo->molecularWeights(); 
// the species mass fractions are the last components in the solution 
// vector, so the total number of components is the number of species 
// plus the offset of the first mass fraction. 
m_nv = c_offset_Y + m_nsp; 
// enable all species equations by default 
m_do_species.resize(m_nsp, true); 
// but turn off the energy equation at all points 
m_do_energy.resize(m_points,false); 
m_diff.resize(m_nsp*m_points); 
m_multidiff.resize(m_nsp*m_nsp*m_points); 
m_flux.resize(m_nsp,m_points); 
m_wdot.resize(m_nsp,m_points, 0.0); 
m_surfdot.resize(m_nsp, 0.0); 
m_ybar.resize(m_nsp); 
//-------------- default solution bounds -------------------- 
vector_fp vmin(m_nv), vmax(m_nv); 
// no bounds on u 
vmin[0] = -1.e20; 
vmax[0] = 1.e20; 
// V 
vmin[1] = -1.e20; 
vmax[1] = 1.e20; 
// temperature bounds 
vmin[2] = 200.0; 
vmax[2]= 1.e9; 
// lamda should be negative 
vmin[3] = -1.e20; 
vmax[3] = 1.e20; 
// mass fraction bounds 
int k; 
for (k = 0; k < m_nsp; k++) { 
vmin[4+k] = -1.0e-5; 
vmax[4+k] = 1.0e5; 
} 
setBounds(vmin.size(), DATA_PTR(vmin), vmax.size(), DATA_PTR(vmax)); 
//-------------------- default error tolerances ---------------- 
vector_fp rtol(m_nv, 1.0e-8); 
vector_fp atol(m_nv, 1.0e-15); 
setTolerances(rtol.size(), DATA_PTR(rtol), atol.size(), DATA_PTR(atol),false); 
setTolerances(rtol.size(), DATA_PTR(rtol), atol.size(), DATA_PTR(atol),true); 
//-------------------- grid refinement ------------------------- 
m_refiner->setActive(0, false); 
m_refiner->setActive(1, false); 
m_refiner->setActive(2, false); 
m_refiner->setActive(3, false); 
vector_fp gr; 
for (int ng = 0; ng < m_points; ng++) gr.push_back(1.0*ng/m_points); 
setupGrid(m_points, DATA_PTR(gr)); 


 setID("stagnation flow"); 
ifstream in("Properties2.txt"); //Read in the solid properties 
double proper; 
in>>proper; 
pore1=proper; 
in>>proper; 
pore2=proper; 
in>>proper; 
diam1=proper; 
in>>proper; 
diam2=proper; 
in>>proper; 
Omega1=proper; 
in>>proper; 
Omega2=proper; 
in>>proper; 
srho=proper; 
in>>proper; 
sCp=proper; 
in.close(); 
} 
/** 
* Change the grid size. Called after grid refinement. 
*/ 
void StFlow::resize(int ncomponents, int points) { 
Domain1D::resize(ncomponents, points); 
m_rho.resize(m_points, 0.0); 
m_wtm.resize(m_points, 0.0); 
m_cp.resize(m_points, 0.0); 
m_enth.resize(m_points, 0.0); 
m_visc.resize(m_points, 0.0); 
m_tcon.resize(m_points, 0.0); 
if (m_transport_option == c_Mixav_Transport) { 
m_diff.resize(m_nsp*m_points); 
} 
else { 
m_multidiff.resize(m_nsp*m_nsp*m_points); 
m_diff.resize(m_nsp*m_points); 
} 
m_flux.resize(m_nsp,m_points); 
m_wdot.resize(m_nsp,m_points, 0.0); 
m_do_energy.resize(m_points,false); 
m_fixedy.resize(m_nsp, m_points); 
m_fixedtemp.resize(m_points); 
m_dz.resize(m_points-1); 
m_z.resize(m_points); 
} 
void StFlow::setupGrid(int n, const doublereal* z) { 
resize(m_nv, n); 
int j; 
m_z[0] = z[0]; 
for (j = 1; j < m_points; j++) { 
m_z[j] = z[j]; 
m_dz[j-1] = m_z[j] - m_z[j-1]; 
} 
} 
/** 
* Install a transport manager. 
*/ 


 void StFlow::setTransport(Transport& trans, bool withSoret) { 
m_trans = &trans; 
m_do_soret = withSoret; 
if (m_trans->model() == cMulticomponent) { 
m_transport_option = c_Multi_Transport; 
m_multidiff.resize(m_nsp*m_nsp*m_points); 
m_diff.resize(m_nsp*m_points); 
m_dthermal.resize(m_nsp, m_points, 0.0); 
} 
else if (m_trans->model() == cMixtureAveraged) { 
m_transport_option = c_Mixav_Transport; 
m_diff.resize(m_nsp*m_points); 
if (withSoret) 
throw CanteraError("setTransport", 
"Thermal diffusion (the Soret effect) " 
"requires using a multicomponent transport model."); 
} 
else 
throw CanteraError("setTransport","unknown transport model."); 
} 
/** 
* Set the gas object state to be consistent with the solution at 
* point j. 
*/ 
void StFlow::setGas(const doublereal* x,int j) { 
m_thermo->setTemperature(T(x,j)); 
const doublereal* yy = x + m_nv*j + c_offset_Y; 
m_thermo->setMassFractions_NoNorm(yy); 
m_thermo->setPressure(m_press); 
} 
/** 
* Set the gas state to be consistent with the solution at the 
* midpoint between j and j + 1. 
*/ 
void StFlow::setGasAtMidpoint(const doublereal* x,int j) { 
m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1))); 
const doublereal* yyj = x + m_nv*j + c_offset_Y; 
const doublereal* yyjp = x + m_nv*(j+1) + c_offset_Y; 
for (int k = 0; k < m_nsp; k++) 
m_ybar[k] = 0.5*(yyj[k] + yyjp[k]); 
m_thermo->setMassFractions_NoNorm(DATA_PTR(m_ybar)); 
m_thermo->setPressure(m_press); 
} 
void StFlow::_finalize(const doublereal* x) { 
int k, j; 
doublereal zz, tt; 
int nz = m_zfix.size(); 
bool e = m_do_energy[0]; 
for (j = 0; j < m_points; j++) { 
if (e || nz == 0) 
setTemperature(j, T(x, j)); 
else { 
zz = (z(j) - z(0))/(z(m_points - 1) - z(0)); 
tt = linearInterp(zz, m_zfix, m_tfix); 
setTemperature(j, tt); 
} 
for (k = 0; k < m_nsp; k++) { 
setMassFraction(j, k, Y(x, k, j)); 
} 
} 
if (e) solveEnergyEqn(); 
} 



 //------------------------------------------------------ 
/** 
* Evaluate the residual function for axisymmetric stagnation 
* flow. If jpt is less than zero, the residual function is 
* evaluated at all grid points. If jpt >= 0, then the residual 
* function is only evaluated at grid points jpt-1, jpt, and 
* jpt+1. This option is used to efficiently evaluate the 
* Jacobian numerically. 
* 
*/ 
void AxiStagnFlow::eval(int jg, doublereal* xg, 
doublereal* rg, integer* diagg, doublereal rdt) { 
// if evaluating a Jacobian, and the global point is outside 
// the domain of influence for this domain, then skip 
// evaluating the residual 
if (jg >=0 && (jg < firstPoint() - 1 || jg > lastPoint() + 1)) return; 
// if evaluating a Jacobian, compute the steady-state residual 
if (jg >= 0) rdt = 0.0; 
// start of local part of global arrays 
doublereal* x = xg + loc(); 
doublereal* rsd = rg + loc(); 
integer* diag = diagg + loc(); 
int jmin, jmax, jpt; 
jpt = jg - firstPoint(); 
if (jg < 0) { // evaluate all points 
jmin = 0; 
jmax = m_points - 1; 
} 
else { // evaluate points for Jacobian 
jmin = max(jpt-1, 0); 
jmax = min(jpt+1,m_points-1); 
} 
// properties are computed for grid points from j0 to j1 
int j0 = max(jmin-1,0); 
int j1 = min(jmax+1,m_points-1); 
int j, k; 
//----------------------------------------------------- 
// update properties 
//----------------------------------------------------- 
// update thermodynamic properties only if a Jacobian is not 
// being evaluated 
if (jpt < 0) { //if (jpt < 0 || (m_transport_option == c_Multi_Transport)) { 
updateThermo(x, j0, j1); 
// update transport properties only if a Jacobian is not being 
// evaluated 
updateTransport(x, j0, j1); 
} 
// update the species diffusive mass fluxes whether or not a 
// Jacobian is being evaluated 
updateDiffFluxes(x, j0, j1); 
//---------------------------------------------------- 
// evaluate the residual equations at all required 
// grid points 
//---------------------------------------------------- 


 doublereal sum, sum2, dtdzj; 
doublereal lam, visc, Re; //Defining new variables. 
double length=m_points; // 
hconv.resize(length); // 
//initialize property vectors 
vector<double> pore(length); 
vector<double> diam(length); 
vector<double> scond(length); 
vector<double> Omega(length); 
vector<double> Cmult(length); 
vector<double> mpow(length); 
vector<double> RK(length); 
//populate property vectors 
for (int i=0; i<=length-1;i++) 
{ 
if (z(i)<0.033) 
{ 
pore[i]=pore1; 
diam[i]=diam1; 
} 
else if (z(i)>0.037) 
{ 
pore[i]=pore2; 
diam[i]=diam2; 
} 
else 
{ 
pore[i]=(((pore2-pore1)/(.037-.033))*(z(i)-0.033))+pore1; 
diam[i]=(((diam2-diam1)/(.037-.033))*(z(i)-0.033))+diam1; 
} 
RK[i]=(3*(1-pore[i])/diam[i]); 
Cmult[i]=-400*diam[i]+0.687; 
mpow[i]=443.7*diam[i]+0.361; 
scond[i]=0.188-17.5*diam[i]; 
} 
for (int i=0; i<=length-1;i++) 
{ 
if (z(i)<0.035) 
{ 
Omega[i]=Omega1; 
} 
else 
{ 
Omega[i]=Omega2; 
} 
} 
int solidenergy=0; 
//loop over gas energy vecotr. If it is going to be solved then find hv 
for(j=jmin;j<=jmax;j++) 
{ 
solidenergy+=m_do_energy[j]; 
} 
solidenergy=1; 
if (solidenergy!=0) 
{ 
for (j = jmin; j <= jmax; j++) 
{ 
lam=m_tcon[j]; 
visc=m_visc[j]; 
Re=(rho_u(x,j)*pore[j]*diam[j])/visc; 
hconv[j]=((lam*Cmult[j]*pow(Re,mpow[j]))/pow(diam[j],2)); 
} 
//Solve for the solid profile if required 
if (dosolid==1) 
{ 


 solid(x,hconv,scond,RK,Omega,srho,sCp,rdt); 
dosolid=0; 
} 
} 
for (j = jmin; j <= jmax; j++) { 
//---------------------------------------------- 
// left boundary 
//---------------------------------------------- 
if (j == 0) { 
// these may be modified by a boundary object 
// Continuity. This propagates information right-to-left, 
// since rho_u at point 0 is dependent on rho_u at point 1, 
// but not on mdot from the inlet. 
rsd[index(c_offset_U,0)] = 
-(rho_u(x,1) - rho_u(x,0))/m_dz[0] 
-(density(1)*V(x,1) + density(0)*V(x,0)); 
// the inlet (or other) object connected to this one 
// will modify these equations by subtracting its values 
// for V, T, and mdot. As a result, these residual equations 
// will force the solution variables to the values for 
// the boundary object 
rsd[index(c_offset_V,0)] = V(x,0); 
rsd[index(c_offset_T,0)] = T(x,0); 
rsd[index(c_offset_L,0)] = -rho_u(x,0); 
// The default boundary condition for species is zero 
// flux. However, the boundary object may modify 
// this. 
sum = 0.0; 
for (k = 0; k < m_nsp; k++) { 
sum += Y(x,k,0); 
rsd[index(c_offset_Y + k, 0)] = 
-(m_flux(k,0) + rho_u(x,0)* Y(x,k,0)); 
} 
rsd[index(c_offset_Y, 0)] = 1.0 - sum; 
} 
//---------------------------------------------- 
// 
// right boundary 
// 
//---------------------------------------------- 
else if (j == m_points - 1) { 
// the boundary object connected to the right of this 
// one may modify or replace these equations. The 
// default boundary conditions are zero u, V, and T, 
// and zero diffusive flux for all species. 
rsd[index(0,j)] = rho_u(x,j); 
rsd[index(1,j)] = V(x,j); 
rsd[index(2,j)] = T(x,j); 
rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1); 
diag[index(c_offset_L, j)] = 0; 
doublereal sum = 0.0; 
for (k = 0; k < m_nsp; k++) { 
sum += Y(x,k,j); 
rsd[index(k+4,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j); 
} 
rsd[index(4,j)] = 1.0 - sum; 
diag[index(4,j)] = 0; 



 } 
//------------------------------------------ 
// interior points 
//------------------------------------------ 
else { 
//---------------------------------------------- 
// Continuity equation 
// 
// Note that this propagates the mass flow rate 
// information to the left (j+1 -> j) from the 
// value specified at the right boundary. The 
// lambda information propagates in the opposite 
// direction. 
// 
// d(\rho u)/dz + 2\rho V = 0 
// 
//------------------------------------------------ 
rsd[index(c_offset_U,j)] = 
-(rho_u(x,j+1)*pore[j+1] - rho_u(x,j)*pore[j])/m_dz[j] //added porosity 
-(density(j+1)*V(x,j+1) + density(j)*V(x,j)); 
//algebraic constraint 
diag[index(c_offset_U, j)] = 0; 
//------------------------------------------------ 
// Radial momentum equation 
// 
// \rho u dV/dz + \rho V^2 = d(\mu dV/dz)/dz - lambda 
// 
//------------------------------------------------- 
rsd[index(c_offset_V,j)] 
= (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j) 
- m_rho[j]*V(x,j)*V(x,j))/m_rho[j] 
- rdt*(V(x,j) - V_prev(j)); 
diag[index(c_offset_V, j)] = 1; 
//------------------------------------------------- 
// Species equations 
// 
// \rho u dY_k/dz + dJ_k/dz + M_k\omega_k 
// 
//------------------------------------------------- 
getWdot(x,j); 
doublereal convec, diffus; 
for (k = 0; k < m_nsp; k++) { 
convec = rho_u(x,j)*dYdz(x,k,j)*pore[j]; //added porosity 
diffus = 2.0*(m_flux(k,j)*pore[j] - m_flux(k,j-1)*pore[j-1]) //added 
porosity 
/(z(j+1) - z(j-1)); 
rsd[index(c_offset_Y + k, j)] 
= (m_wt[k]*(wdot(k,j)*pore[j] ) //added porosity 
- convec - diffus)/(m_rho[j]*pore[j]) //added porosity 
- rdt*(Y(x,k,j) - Y_prev(k,j)); 
diag[index(c_offset_Y + k, j)] = 1; 
} 
//----------------------------------------------- 
// energy equation 
//----------------------------------------------- 
if (m_do_energy[j]) { 



 setGas(x,j); 
// heat release term 
const vector_fp& h_RT = m_thermo->enthalpy_RT_ref(); 
const vector_fp& cp_R = m_thermo->cp_R_ref(); 
sum = 0.0; 
sum2 = 0.0; 
doublereal flxk; 
for (k = 0; k < m_nsp; k++) { 
flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j)); 
sum += wdot(k,j)*h_RT[k]; 
sum2 += flxk*cp_R[k]/m_wt[k]; 
} 
sum *= GasConstant * T(x,j); 
dtdzj = dTdz(x,j); 
sum2 *= GasConstant * dtdzj; 
rsd[index(c_offset_T, j)] = 
- m_cp[j]*rho_u(x,j)*dtdzj 
- divHeatFlux(x,j) - sum - sum2; 
rsd[index(c_offset_T, j)] = 
//adding of convective term 
rsd[index(c_offset_T, j)] - (hconv[j]*(T(x,j)-
Tw[j]))/pore[j]; // 
rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]); 
rsd[index(c_offset_T, j)] = 
rsd[index(c_offset_T, j)] + m_efctr*(T_fixed(j) - T(x,j)); 
rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j)); 
diag[index(c_offset_T, j)] = 1; 
} 
// residual equations if the energy equation is disabled 
if (!m_do_energy[j]) { 
rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j); 
diag[index(c_offset_T, j)] = 0; 
} 
rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1); 
diag[index(c_offset_L, j)] = 0; 
} 
} 
} 
//Solid solver 
void StFlow::solid(doublereal* x, vector<double> &hconv, vector<double>& scond, 
vector<double>& RK, vector<double>&Omega,double & srho,double & sCp, double rdt) { 
writelog("Computing Solid Temperature Field..."); 
double length=m_points; // 
Tw.resize(length); 
Twprev.resize(length); 
dq.resize(length); 
if (adapt==1) 
{ 
int j=0; 
//Ensures solid profile has the correct number of points after an adaption 
for (int i=0;i<=length-1;i++) 
{ 
if (z(i)==zprev[j]) //same point 
{ 
Twprev[i]=Twprev1[j]; 
j++; 
} 
else if (z(i)>zprev[j]) //deleted point 


 { 
i--; 
j++; 
} 
else if (z(i)<zprev[j]) //added point 
{ 
if (i==1) 
{ 
Twprev[i]=Twprev[i-1]; 
} 
else 
{ 
Twprev[i]=Twprev[i-1]+(Twprev[i-1]-Twprev[i-2])*((z(i)-z(i-
1))/(z(i-1)-z(i-2))); 
} 
} 
} 
} 
adapt=1; 
zprev.resize(length); 
Twprev1.resize(length); 
//Start of Conduction Radiation Stuff 
// 
//Vector Iinitialization 
vector<double> edia(length); 
vector<double> fdia(length); 
vector<double> gdia(length); 
vector<double> rhs(length); 
vector<double> dqnew(length); 
double sigma=0.0000000567; 
double change1=1; 
//Vector Population 
for(int i=0;i<=length-1;i++) 
{ 
dq[i]=0; 
} 
double T0=300; 
double T1=300; 
int count1=0; 
int fail1=0; 
while (change1>0.000001) 
{ 
count1=count1+1; 
for(int i=0;i<=length-1;i++) 
{ 
if (i==0) 
{ 
edia[i]=0; 
fdia[i]=1; 
gdia[i]=-1; 
rhs[i]=0; 
} 
else if (i==length-1) 
{ 
edia[i]=-1; 
fdia[i]=1; 
gdia[i]=0; 
rhs[i]=0; 
} 
else 
{ 
edia[i]=(2*scond[i])/((z(i)-z(i-1))*(z(i+1)-z(i-1))); 
fdia[i]=-(2*scond[i])/((z(i+1)-z(i))*(z(i+1)-z(i-1)))-
(2*scond[i])/((z(i)-z(i-1))*(z(i+1)-z(i-1)))-hconv[i]-srho*sCp*rdt; 
gdia[i]=(2*scond[i])/((z(i+1)-z(i))*(z(i+1)-z(i-1))); 
rhs[i]=-hconv[i]*T(x,i)+dq[i]-srho*sCp*rdt*Twprev[i]; 
} 
} 



 //Decomposition 
for(int i=1;i<=length-1;i++) 
{ 
edia[i]=edia[i]/fdia[i-1]; 
fdia[i]=fdia[i]-edia[i]*gdia[i-1]; 
} 
//Forward Substitution 
for(int i=1;i<=length-1;i++) 
{ 
rhs[i]=rhs[i]-edia[i]*rhs[i-1]; 
} 
//Back Substitution 
Tw[length-1]=rhs[length-1]/fdia[length-1]; 
for(int i=length-2;i>=0;i--) 
{ 
Tw[i]=(rhs[i]-gdia[i]*Tw[i+1])/fdia[i]; 
} 
T0=Tw[0]; 
T1=Tw[length-1]; 
//Radiation Time 
//Vector Initialization 
vector<double> qplus(length); 
vector<double> qpnew(length); 
vector<double> qminus(length); 
vector<double> qmnew(length); 
double change2=1; 
//Vector Population 
for(int i=0;i<=length-1;i++) 
{ 
double temp=T(x,i); 
double temp2; 
if (i==0) 
{ 
qplus[i]=sigma*pow(temp,4); 
qpnew[i]=sigma*pow(temp,4); 
qminus[i]=0; 
qmnew[i]=0; 
temp2=temp; 
} 
else if (i==length-1) 
{ 
qplus[i]=0; 
qpnew[i]=0; 
qminus[i]=sigma*pow(temp2,4); 
qmnew[i]=sigma*pow(temp2,4); 
} 
else 
{ 
qplus[i]=0; 
qpnew[i]=0; 
qminus[i]=0; 
qmnew[i]=0; 
} 
} 
int count=0; 
int fail=0; 
S2 method 
while (change2>0.000001) 
{ 
count=count+1; 
for(int i=1;i<=length-1;i++) 
{ 
double temp=Tw[i]; 
qpnew[i]=(qpnew[i-1]+RK[i]*(z(i)-z(i-
1))*Omega[i]*qminus[i]+2*RK[i]*(z(i)-z(i-1))*(1-Omega[i])*sigma*pow(temp,4))/(1+(z(i)-z(i-
1))*RK[i]*(2-Omega[i])); 
} 
for(int i=length-2;i>=0;i--) 


 { 
double temp=Tw[i]; 
qmnew[i]=(qmnew[i+1]+RK[i]*(z(i+1)-
z(i))*Omega[i]*qpnew[i]+2*RK[i]*(z(i+1)-z(i))*(1-Omega[i])*sigma*pow(temp,4))/(1+(z(i+1)-
z(i))*RK[i]*(2-Omega[i])); 
} 
double norm1=0; 
double norm2=0; 
for(int i=0;i<=length-1;i++) 
{ 
norm1+=(qpnew[i]-qplus[i])*(qpnew[i]-qplus[i]); 
norm2+=(qmnew[i]-qminus[i])*(qmnew[i]-qminus[i]); 
qplus[i]=qpnew[i]; 
qminus[i]=qmnew[i]; 
} 
norm1=sqrt(norm1); 
norm2=sqrt(norm2); 
if (count>100) 
{ 
change2=0; 
fail=1; 
} 
else 
{ 
change2=max(norm1,norm2); 
} 
} 
if (fail==1) 
{ 
for(int i=0;i<=length-1;i++) 
{ 
dqnew[i]=dq[i]; 
} 
writelog("Rad Stall"); 
} 
else 
{ 
for(int i=0;i<=length-1;i++) 
{ 
double temp=Tw[i]; 
dqnew[i]=4*RK[i]*(1-Omega[i])*(sigma*pow(temp,4)-
0.5*qplus[i]-0.5*qminus[i]); 
} 
} 
double norm=0; 
double a=0.1; 
for (int i=0;i<=length-1;i++) 
{ 
norm+=(dqnew[i]-dq[i])*(dqnew[i]-dq[i]); 
dq[i]=a*dqnew[i]+(1-a)*dq[i]; 
} 
if (count1>400) 
{ 
fail1=1; 
change1=0; 
} 
else 
{ 
change1=sqrt(norm); 
} 
} 
if (fail1==1) 
{ 
for (int i=0;i<=length-1;i++) 
{ 
Tw[i]=Twprev1[i]; 
} 
writelog("Rad not Converged"); 
} 
for (int i=0;i<=length-1;i++) 
{ 


 Twprev1[i]=Tw[i]; 
zprev[i]=z(i); 
} 
writelog("Success\n"); 
// 
//End of Newly Added Conduction Radiation Stuff 
} 
/** 
* Update the transport properties at grid points in the range 
* from j0 to j1, based on solution x. 
*/ 
void StFlow::updateTransport(doublereal* x,int j0, int j1) { 
int j,k,m; 
if (m_transport_option == c_Mixav_Transport) { 
for (j = j0; j < j1; j++) { 
setGasAtMidpoint(x,j); 
m_visc[j] = m_trans->viscosity(); 
m_trans->getMixDiffCoeffs(DATA_PTR(m_diff) + j*m_nsp); 
m_tcon[j] = m_trans->thermalConductivity(); 
} 
} 
else if (m_transport_option == c_Multi_Transport) { 
doublereal sum, sumx, wtm, dz; 
doublereal eps = 1.0e-12; 
for (m = j0; m < j1; m++) { 
setGasAtMidpoint(x,m); 
dz = m_z[m+1] - m_z[m]; 
wtm = m_thermo->meanMolecularWeight(); 
m_visc[j] = m_trans->viscosity(); 
m_trans->getMultiDiffCoeffs(m_nsp, 
DATA_PTR(m_multidiff) + mindex(0,0,m)); 
for (k = 0; k < m_nsp; k++) { 
sum = 0.0; 
sumx = 0.0; 
for (j = 0; j < m_nsp; j++) { 
if (j != k) { 
sum += m_wt[j]*m_multidiff[mindex(k,j,m)]* 
((X(x,j,m+1) - X(x,j,m))/dz + eps); 
sumx += (X(x,j,m+1) - X(x,j,m))/dz; 
} 
} 
m_diff[k + m*m_nsp] = sum/(wtm*(sumx+eps)); 
} 
m_tcon[m] = m_trans->thermalConductivity(); 
if (m_do_soret) { 
m_trans->getThermalDiffCoeffs(m_dthermal.ptrColumn(0) + m*m_nsp); 
} 
} 
} 
} 
//------------------------------------------------------ 
/** 
* Evaluate the residual function for axisymmetric stagnation 
* flow. If jpt is less than zero, the residual function is 
* evaluated at all grid points. If jpt >= 0, then the residual 
* function is only evaluated at grid points jpt-1, jpt, and 
* jpt+1. This option is used to efficiently evaluate the 
* Jacobian numerically. 
* 
*/ 



 void FreeFlame::eval(int jg, doublereal* xg, 
doublereal* rg, integer* diagg, doublereal rdt) { 
// if evaluating a Jacobian, and the global point is outside 
// the domain of influence for this domain, then skip 
// evaluating the residual 
if (jg >=0 && (jg < firstPoint() - 1 || jg > lastPoint() + 1)) return; 
// if evaluating a Jacobian, compute the steady-state residual 
if (jg >= 0) rdt = 0.0; 
// start of local part of global arrays 
doublereal* x = xg + loc(); 
doublereal* rsd = rg + loc(); 
integer* diag = diagg + loc(); 
int jmin, jmax, jpt; 
jpt = jg - firstPoint(); 
if (jg < 0) { // evaluate all points 
jmin = 0; 
jmax = m_points - 1; 
} 
else { // evaluate points for Jacobian 
jmin = max(jpt-1, 0); 
jmax = min(jpt+1,m_points-1); 
} 
// properties are computed for grid points from j0 to j1 
int j0 = max(jmin-1,0); 
int j1 = min(jmax+1,m_points-1); 
int j, k; 
//----------------------------------------------------- 
// update properties 
//----------------------------------------------------- 
// update thermodynamic properties only if a Jacobian is not 
// being evaluated 
if (jpt < 0) { 
updateThermo(x, j0, j1); 
updateTransport(x, j0, j1); 
} 
// update the species diffusive mass fluxes whether or not a 
// Jacobian is being evaluated 
updateDiffFluxes(x, j0, j1); 
//---------------------------------------------------- 
// evaluate the residual equations at all required 
// grid points 
//---------------------------------------------------- 
doublereal sum, sum2, dtdzj; 
for (j = jmin; j <= jmax; j++) { 
//---------------------------------------------- 
// left boundary 
//---------------------------------------------- 
if (j == 0) { 
// these may be modified by a boundary object 
// Continuity. This propagates information right-to-left, 


 // since rho_u at point 0 is dependent on rho_u at point 1, 
// but not on mdot from the inlet. 
rsd[index(c_offset_U,0)] = 
-(rho_u(x,1) - rho_u(x,0))/m_dz[0] 
-(density(1)*V(x,1) + density(0)*V(x,0)); 
// the inlet (or other) object connected to this one 
// will modify these equations by subtracting its values 
// for V, T, and mdot. As a result, these residual equations 
// will force the solution variables to the values for 
// the boundary object 
rsd[index(c_offset_V,0)] = V(x,0); 
rsd[index(c_offset_T,0)] = T(x,0); 
rsd[index(c_offset_L,0)] = -rho_u(x,0); 
// The default boundary condition for species is zero 
// flux 
sum = 0.0; 
for (k = 0; k < m_nsp; k++) { 
sum += Y(x,k,0); 
rsd[index(c_offset_Y + k, 0)] = 
-(m_flux(k,0) + rho_u(x,0)* Y(x,k,0)); 
} 
rsd[index(c_offset_Y, 0)] = 1.0 - sum; 
} 
//---------------------------------------------- 
// 
// right boundary 
// 
//---------------------------------------------- 
else if (j == m_points - 1) { 
// the boundary object connected to the right of this 
// one may modify or replace these equations. The 
// default boundary conditions are zero u, V, and T, 
// and zero diffusive flux for all species. 
// zero gradient 
rsd[index(0,j)] = rho_u(x,j) - rho_u(x,j-1); 
rsd[index(1,j)] = V(x,j); 
rsd[index(2,j)] = T(x,j) - T(x,j-1); 
doublereal sum = 0.0; 
rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1); 
diag[index(c_offset_L, j)] = 0; 
for (k = 0; k < m_nsp; k++) { 
sum += Y(x,k,j); 
rsd[index(k+4,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j); 
} 
rsd[index(4,j)] = 1.0 - sum; 
diag[index(4,j)] = 0; 
} 
//------------------------------------------ 
// interior points 
//------------------------------------------ 
else { 
//---------------------------------------------- 
// Continuity equation 
//---------------------------------------------- 
if (grid(j) > m_zfixed){ 
rsd[index(c_offset_U,j)] = 
- (rho_u(x,j) - rho_u(x,j-1))/m_dz[j-1] 
- (density(j-1)*V(x,j-1) + density(j)*V(x,j)); 
} 



 else if (grid(j) == m_zfixed){ 
if (m_do_energy[j]) { 
rsd[index(c_offset_U,j)] = (T(x,j) - m_tfixed); 
} 
else { 
rsd[index(c_offset_U,j)] = (rho_u(x,j) 
- m_rho[0]*0.3); 
} 
} 
else if(grid(j) < m_zfixed){ 
rsd[index(c_offset_U,j)] = 
- (rho_u(x,j+1) - rho_u(x,j))/m_dz[j] 
- (density(j+1)*V(x,j+1) + density(j)*V(x,j)); 
} 
//algebraic constraint 
diag[index(c_offset_U, j)] = 0; 
//------------------------------------------------ 
// Radial momentum equation 
// 
// \rho u dV/dz + \rho V^2 = d(\mu dV/dz)/dz - lambda 
// 
//------------------------------------------------- 
rsd[index(c_offset_V,j)] 
= (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j) 
- m_rho[j]*V(x,j)*V(x,j))/m_rho[j] 
- rdt*(V(x,j) - V_prev(j)); 
diag[index(c_offset_V, j)] = 1; 
//------------------------------------------------- 
// Species equations 
// 
// \rho u dY_k/dz + dJ_k/dz + M_k\omega_k 
// 
//------------------------------------------------- 
getWdot(x,j); 
doublereal convec, diffus; 
for (k = 0; k < m_nsp; k++) { 
convec = rho_u(x,j)*dYdz(x,k,j); 
diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1)) 
/(z(j+1) - z(j-1)); 
rsd[index(c_offset_Y + k, j)] 
= (m_wt[k]*(wdot(k,j) ) 
- convec - diffus)/m_rho[j] 
- rdt*(Y(x,k,j) - Y_prev(k,j)); 
diag[index(c_offset_Y + k, j)] = 1; 
} 
//----------------------------------------------- 
// energy equation 
//----------------------------------------------- 
if (m_do_energy[j]) { 
setGas(x,j); 
// heat release term 
const vector_fp& h_RT = m_thermo->enthalpy_RT_ref(); 
const vector_fp& cp_R = m_thermo->cp_R_ref(); 
sum = 0.0; 
sum2 = 0.0; 
doublereal flxk; 
for (k = 0; k < m_nsp; k++) { 
flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j)); 
sum += wdot(k,j)*h_RT[k]; 
sum2 += flxk*cp_R[k]/m_wt[k]; 
} 


 sum *= GasConstant * T(x,j); 
dtdzj = dTdz(x,j); 
sum2 *= GasConstant * dtdzj; 
rsd[index(c_offset_T, j)] = 
- m_cp[j]*rho_u(x,j)*dtdzj 
- divHeatFlux(x,j) - sum - sum2; 
rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]); 
rsd[index(c_offset_T, j)] = 
rsd[index(c_offset_T, j)] + m_efctr*(T_fixed(j) - T(x,j)); 
rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j)); 
diag[index(c_offset_T, j)] = 1; 
} 
// residual equations if the energy equation is disabled 
else { 
rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j); 
diag[index(c_offset_T, j)] = 0; 
} 
rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1); 
diag[index(c_offset_L, j)] = 0; 
} 
} 
} 
/** 
* Print the solution. 
*/ 
void StFlow::showSolution(const doublereal* x) { 
int nn = m_nv/5; 
int i, j, n; 
//char* buf = new char[100]; 
char buf[100]; 
// The mean molecular weight is needed to convert 
updateThermo(x, 0, m_points-1); 
sprintf(buf, " Pressure: %10.4g Pa \n", m_press); 
writelog(buf); 
for (i = 0; i < nn; i++) { 
st_drawline(); 
sprintf(buf, "\n z "); 
writelog(buf); 
for (n = 0; n < 5; n++) { 
sprintf(buf, " %10s ",componentName(i*5 + n).c_str()); 
writelog(buf); 
} 
st_drawline(); 
for (j = 0; j < m_points; j++) { 
sprintf(buf, "\n %10.4g ",m_z[j]); 
writelog(buf); 
for (n = 0; n < 5; n++) { 
sprintf(buf, " %10.4g ",component(x, i*5+n,j)); 
writelog(buf); 
} 
} 
writelog("\n"); 
} 
int nrem = m_nv - 5*nn; 
st_drawline(); 
sprintf(buf, "\n z "); 
writelog(buf); 
for (n = 0; n < nrem; n++) { 
sprintf(buf, " %10s ", componentName(nn*5 + n).c_str()); 
writelog(buf); 
} 
st_drawline(); 
for (j = 0; j < m_points; j++) { 


 sprintf(buf, "\n %10.4g ",m_z[j]); 
writelog(buf); 
for (n = 0; n < nrem; n++) { 
sprintf(buf, " %10.4g ",component(x, nn*5+n,j)); 
writelog(buf); 
} 
} 
writelog("\n"); 
} 
/** 
* Update the diffusive mass fluxes. 
*/ 
void StFlow::updateDiffFluxes(const doublereal* x, int j0, int j1) { 
int j, k, m; 
doublereal sum, wtm, rho, dz, gradlogT; 
switch (m_transport_option) { 
case c_Mixav_Transport: 
case c_Multi_Transport: 
for (j = j0; j < j1; j++) { 
sum = 0.0; 
wtm = m_wtm[j]; 
rho = density(j); 
dz = z(j+1) - z(j); 
for (k = 0; k < m_nsp; k++) { 
m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm); 
m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz; 
sum -= m_flux(k,j); 
} 
// correction flux to insure that \sum_k Y_k V_k = 0. 
for (k = 0; k < m_nsp; k++) m_flux(k,j) += sum*Y(x,k,j); 
} 
break; 
default: 
throw CanteraError("updateDiffFluxes","unknown transport model"); 
} 
if (m_do_soret) { 
for (m = j0; m < j1; m++) { 
gradlogT = 2.0*(T(x,m+1) - T(x,m))/(T(x,m+1) + T(x,m)); 
for (k = 0; k < m_nsp; k++) { 
m_flux(k,m) -= m_dthermal(k,m)*gradlogT; 
} 
} 
} 
} 
string StFlow::componentName(int n) const { 
switch(n) { 
case 0: return "u"; 
case 1: return "V"; 
case 2: return "T"; 
case 3: return "lambda"; 
default: 
if (n >= (int) c_offset_Y && n < (int) (c_offset_Y + m_nsp)) { 
return m_thermo->speciesName(n - c_offset_Y); 
} 
else 
return "<unknown>"; 
} 
} 
//added by Karl Meredith 
int StFlow::componentIndex(string name) const { 


 if(name=="u") {return 0;} 
else if (name=="V") {return 1;} 
else if (name=="T") {return 2;} 
else if (name=="lambda") {return 3;} 
else { 
for (int n=4;n<m_nsp+4;n++){ 
if(componentName(n)==name){ 
return n; 
} 
} 
} 
return -1; 
} 
void StFlow::restore(const XML_Node& dom, doublereal* soln) { 
vector<string> ignored; 
int nsp = m_thermo->nSpecies(); 
vector_int did_species(nsp, 0); 
vector<XML_Node*> str; 
dom.getChildren("string",str); 
int nstr = static_cast<int>(str.size()); 
for (int istr = 0; istr < nstr; istr++) { 
const XML_Node& nd = *str[istr]; 
writelog(nd["title"]+": "+nd.value()+"\n"); 
} 
//map<string, double> params; 
double pp = -1.0; 
pp = getFloat(dom, "pressure", "pressure"); 
setPressure(pp); 
vector<XML_Node*> d; 
dom.child("grid_data").getChildren("floatArray",d); 
int nd = static_cast<int>(d.size()); 
vector_fp x; 
int n, np = 0, j, ks, k; 
string nm; 
bool readgrid = false, wrote_header = false; 
for (n = 0; n < nd; n++) { 
const XML_Node& fa = *d[n]; 
nm = fa["title"]; 
if (nm == "z") { 
getFloatArray(fa,x,false); 
np = x.size(); 
writelog("Grid contains "+int2str(np)+ 
" points.\n"); 
readgrid = true; 
setupGrid(np, DATA_PTR(x)); 
} 
} 
if (!readgrid) { 
throw CanteraError("StFlow::restore", 
"domain contains no grid points."); 
} 
writelog("Importing datasets:\n"); 
for (n = 0; n < nd; n++) { 
const XML_Node& fa = *d[n]; 
nm = fa["title"]; 
getFloatArray(fa,x,false); 
if (nm == "u") { 
writelog("axial velocity "); 
if ((int) x.size() == np) { 


 for (j = 0; j < np; j++) { 
soln[index(0,j)] = x[j]; 
} 
} 
else { 
goto error; 
} 
} 
else if (nm == "z") { 
; // already read grid 
} 
else if (nm == "V") { 
writelog("radial velocity "); 
if ((int) x.size() == np) { 
for (j = 0; j < np; j++) 
soln[index(1,j)] = x[j]; 
} 
else goto error; 
} 
else if (nm == "T") { 
writelog("temperature "); 
if ((int) x.size() == np) { 
for (j = 0; j < np; j++) 
soln[index(2,j)] = x[j]; 
// For fixed-temperature simulations, use the 
// imported temperature profile by default. If 
// this is not desired, call setFixedTempProfile 
// *after* restoring the solution. 
vector_fp zz(np); 
for (int jj = 0; jj < np; jj++) 
zz[jj] = (grid(jj) - zmin())/(zmax() - zmin()); 
setFixedTempProfile(zz, x); 
} 
else goto error; 
} 
else if (nm == "L") { 
writelog("lambda "); 
if ((int) x.size() == np) { 
for (j = 0; j < np; j++) 
soln[index(3,j)] = x[j]; 
} 
else goto error; 
} 
else if (m_thermo->speciesIndex(nm) >= 0) { 
writelog(nm+" "); 
if ((int) x.size() == np) { 
k = m_thermo->speciesIndex(nm); 
did_species[k] = 1; 
for (j = 0; j < np; j++) 
soln[index(k+4,j)] = x[j]; 
} 
} 
else 
ignored.push_back(nm); 
} 
if (ignored.size() != 0) { 
writelog("\n\n"); 
writelog("Ignoring datasets:\n"); 
int nn = static_cast<int>(ignored.size()); 
for (int n = 0; n < nn; n++) { 
writelog(ignored[n]+" "); 
} 
} 
for (ks = 0; ks < nsp; ks++) { 
if (did_species[ks] == 0) { 
if (!wrote_header) { 
writelog("Missing data for species:\n"); 


 wrote_header = true; 
} 
writelog(m_thermo->speciesName(ks)+" "); 
} 
} 
return; 
error: 
throw CanteraError("StFlow::restore","Data size error"); 
} 
void StFlow::save(XML_Node& o, doublereal* sol) { 
int k; 
ArrayViewer soln(m_nv, m_points, sol + loc()); 
//Added to provide Tw as an output 
ofstream f("Solid.txt"); 
ofstream f2("Rad.txt"); 
ofstream f3("Tout.txt"); 
double length=Tw.size(); 
for (int i=0;i<=length-1;i++) 
{ 
f<<Tw[i]<<endl; 
f2<<dq[i]<<endl; 
} 
f3<<Tw[length-1]<<endl; 
f.close(); 
f2.close(); 
f3.close(); 
// 
XML_Node& flow = (XML_Node&)o.addChild("domain"); 
flow.addAttribute("type",flowType()); 
flow.addAttribute("id",m_id); 
flow.addAttribute("points",m_points); 
flow.addAttribute("components",m_nv); 
if (m_desc != "") addString(flow,"description",m_desc); 
XML_Node& gv = flow.addChild("grid_data"); 
addFloat(flow, "pressure", m_press, "Pa", "pressure"); 
addFloatArray(gv,"z",m_z.size(),DATA_PTR(m_z), 
"m","length"); 
vector_fp x(static_cast<size_t>(soln.nColumns())); 
soln.getRow(0,DATA_PTR(x)); 
addFloatArray(gv,"u",x.size(),DATA_PTR(x),"m/s","velocity"); 
soln.getRow(1,DATA_PTR(x)); 
addFloatArray(gv,"V", 
x.size(),DATA_PTR(x),"1/s","rate"); 
soln.getRow(2,DATA_PTR(x)); 
addFloatArray(gv,"T",x.size(),DATA_PTR(x),"K","temperature",0.0); 
soln.getRow(3,DATA_PTR(x)); 
addFloatArray(gv,"L",x.size(),DATA_PTR(x),"N/m^4"); 
for (k = 0; k < m_nsp; k++) { 
soln.getRow(4+k,DATA_PTR(x)); 
addFloatArray(gv,m_thermo->speciesName(k), 
x.size(),DATA_PTR(x),"","massFraction",0.0,1.0); 
} 
} 
void StFlow::setJac(MultiJac* jac) { 
m_jac = jac; 
} 


} // namespace