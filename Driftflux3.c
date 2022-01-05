/* user defined function implementing the drift-flux model*/
// pay attention to the sellting velocity direction in y or z direction

#include "udf.h"
#define dp3 30.0e-7                                     /* need modification */
#define rhoparticle 1000.0
#define miu 1.7894e-5
#define rho 1.225
#define vk (miu/rho)
#define Cc3 (1.0+(0.066e-6/dp3)*(2.34+1.05*exp(-0.39*dp3/0.066e-6)))
#define Brownian3 (1.38e-23*298.0*Cc3/(3.0*3.1415926*miu*dp3))
#define Sc3 (miu/(rho*Brownian3))
#define Sc33 (pow(Sc3,(-1.0/3.0)))
#define vs_x 0.0
#define vs_y 0.0
#define vs_z3 (-rhoparticle*dp3*dp3*9.81*Cc3/(18.0*miu))












DEFINE_UDS_FLUX(particle_flux3,f,t,i)
{
real NV_VEC(vs_vec),NV_VEC(A);
real settling_flux;
real total_flux;
real density;
cell_t c0,c1;
Thread *t0,*t1;

c0=F_C0(f,t);
t0=F_C0_THREAD(f,t);
F_AREA(A,f,t);

NV_D(vs_vec,=,0.0,0.0,vs_z3);                        /* need modification */
settling_flux=NV_DOT(vs_vec,A);

if (BOUNDARY_FACE_THREAD_P(t))
{
if (NNULLP(THREAD_STORAGE(t,SV_DENSITY)))
density=F_R(f,t);
else
density=C_R(c0,t0);
total_flux=F_FLUX(f,t)+settling_flux*density;
}

else
{
c1=F_C1(f,t);
t1=F_C1_THREAD(f,t);
density=(C_R(c0,t0)+C_R(c1,t1))/2.0;
total_flux=F_FLUX(f,t)+settling_flux*density;
}

return total_flux;
}













real deposition_velocity3();
DEFINE_PROFILE(wall_bc3,t,position)
{
face_t f;
cell_t c0;
Thread *t0=t->t0;
real xw[ND_ND];
real xc[ND_ND];
real dx[ND_ND];
real dy;
real cwall;
real vd;
real NV_VEC(vs_vec);
real NV_VEC(A);
real uconst;
real wss;
real prantl;


begin_f_loop(f,t)
{
c0=F_C0(f,t);
F_CENTROID(xw,f,t);
C_CENTROID(xc,c0,t0);
NV_VV(dx, =, xc,-, xw);
dy=ND_MAG(dx[0], dx[1], dx[2]);
prantl=1.0;

NV_D(vs_vec,=,0.0,0.0,vs_z3);                      /* need modification */
F_AREA(A,f,t);
uconst=NV_DOT(vs_vec,A)/NV_MAG(A);
wss=fabs(F_STORAGE_R_N3V(f,t,SV_WALL_SHEAR)[0]+F_STORAGE_R_N3V(f,t,SV_WALL_SHEAR)[1])/NV_MAG(A);
vd=deposition_velocity3(uconst,wss);

if (uconst==0.0)
cwall=C_UDSI(c0,t0,0)-C_UDSI(c0,t0,0)*vd*dy*rho*prantl/C_MU_EFF(c0,t0);
else if (uconst<0.0)
cwall=C_UDSI(c0,t0,0)-C_UDSI(c0,t0,0)*(vd+fabs(vs_z3))*dy*rho*prantl/C_MU_EFF(c0,t0);
else
cwall=C_UDSI(c0,t0,0)-C_UDSI(c0,t0,0)*(vd-fabs(vs_z3))*dy*rho*prantl/C_MU_EFF(c0,t0);

if(cwall>0.0)
F_PROFILE(f,t,position)=cwall;
else
F_PROFILE(f,t,position)=0.0;
}
end_f_loop(f,t)
}


real deposition_velocity3(uconst,wss)
real uconst,wss;
{
real ustar,Rplus,vs;
real A,B,AI,term1,index1,index2,index3,index4,index5,index6,index7,index8;
real vd;
ustar=sqrt(wss/rho);
Rplus=dp3*ustar/(2.0*vk);
vs=uconst/ustar;
index1=pow((10.92*Sc33+4.3),3.0);
index2=1.0/Sc3+0.0609;
index3=8.6-10.92*Sc33;
index4=sqrt(3.0)*10.92*Sc33;
A=0.5*log(index1/index2)+sqrt(3.0)*atan(index3/index4);
index5=pow((10.92*Sc33+Rplus),3.0);
index6=1.0/Sc3+7.669e-4*pow(Rplus,3.0);
index7=2.0*Rplus-10.92*Sc33;
index8=sqrt(3.0)*10.92*Sc33;
B=0.5*log(index5/index6)+sqrt(3.0)*atan(index7/index8);
AI=3.64*pow(Sc3,(2.0/3.0))*(A-B)+39.0;

if(uconst==0.0)
vd=ustar/AI;

else if (uconst>0.0)
{
term1=fabs(vs_z3)*AI/ustar;
if(term1>300.0)
term1=300.0;
vd=fabs(vs_z3)/(1.0-exp(-1.0*term1));
}

else
{
term1=fabs(vs_z3)*AI/ustar;
if(term1>300.0)
term1=300.0;
vd=fabs(vs_z3)/(exp(term1)-1.0);
}

return vd;
}












DEFINE_DIFFUSIVITY(particle_phase_diffusivity3,c,t,i)
{
return (C_MU_EFF(c,t));
}





