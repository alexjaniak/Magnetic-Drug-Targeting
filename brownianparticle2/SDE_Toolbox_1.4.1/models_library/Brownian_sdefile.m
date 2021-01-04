function [out1,out2,out3] = Brownian_sdefile(t,x,flag,bigtheta,SDETYPE,NUMDEPVARS,NUMSIM)
% BROWNIAN_SDEFILE calculates the drift and diffusion vectors of SDE arised from ferrocluster problem
%
%
% [out1,out2,out3] = M9_sdefile(t,x,flag,bigtheta,SDETYPE,NUMDEPVARS,NUMSIM)
%
% IN:     t; working value of independent variable (time)
%         x; working value of dependent variable 
%         flag; a switch, with values 'init' or otherwise
%         bigtheta; complete structural parameter vector
%         SDETYPE; the SDE definition: can be 'Ito' or 'Strat' (Stratonovich)
%         NUMDEPVARS; the number of dependent variables, i.e. the SDE dimension
%         NUMSIM; the number of desired simulations for the SDE numerical integration 
% OUT:    out1; in case of flag='init' is just the initial time, otherwise it is the (vector of) SDE drift(s)
%         out2; in case of flag='init' is the initial value of the dependent variables. Otherwise it is the SDE diffusion(s)
%         out3; in case of flag='init' it is nothing. Otherwise it is the SDE's partial derivative(s) of the diffusion term 
% WORKING:
%		DIM, ii, rparticle, rho, mu0, kB, Temp, chi,rpipe, vmax, vosc,
%		omega,xdiple(1:3),m(1:3), c0, c1, c2, drift, diffusion

%
if(NUMDEPVARS~=3)
    error('Brown_sdefile: NUMDEPVARs is %d while it should be 3', NUMDEPVARS);
end
% Parameters
DIM=NUMDEPVARS;	%dimension of the problem
ii=1;
x0(1:DIM)=bigtheta(ii:ii+DIM-1);	% initial position of cluster
ii=ii+DIM;
xdipole(1:DIM)=bigtheta(ii:ii+DIM-1);	% position of target
ii=ii+DIM;
m(1:DIM)=bigtheta(ii:ii+DIM-1);	% target dipole moment
ii=ii+DIM;
%
rparticle=bigtheta(ii);	% particle radius
ii=ii+1;
rho=bigtheta(ii);	% particle density
ii=ii+1;
mu0=bigtheta(ii);	% vacuum permeability
ii=ii+1;
kB=bigtheta(ii);	% Boltzmann constant
ii=ii+1;
Temp=bigtheta(ii);	% Temperature (Kelvin)
ii=ii+1;
chi=bigtheta(ii);	% suspectibility of particle
ii=ii+1;
vis=bigtheta(ii);	% blood viscosity
ii=ii+1;
rpipe=bigtheta(ii);	% blood vessel radius
ii=ii+1;
vmax=bigtheta(ii);	% blood velocity at vessel center
ii=ii+1;
vosc=bigtheta(ii);	% amplitude of velocity oscillation
ii=ii+1;
omega=bigtheta(ii);	% angular frequency of oscillation
ii=ii+1;
%


if nargin < 3 || isempty(flag)

	c0=6*pi*vis*rparticle;
	c1=3*chi/(3+chi)*mu0/c0*(4*pi/3*rparticle^3);
	c2=sqrt(2.*kB*Temp/c0);

	xparticle=zeros(1,DIM);
	drift=zeros(1,DIM);
	diffusion=zeros(1,DIM);
	out1 = zeros(1,NUMDEPVARS*NUMSIM);
	out2 = zeros(1,NUMDEPVARS*NUMSIM);
	out3 = zeros(1,NUMDEPVARS*NUMSIM);
	for ii=1:NUMSIM;
		ii1=(ii-1)*DIM+1;
		ii2=ii*DIM;
		xparticle(1:DIM)=x(ii1:ii2);
		r=xparticle-xdipole;
%
		velblood=zeros(1,3);
		velblood(1)=poiseuille(xparticle,rpipe,vmax,vosc,omega,t);
		H0=dipoleH(m,r);
		dH=graddipH(m,r);
		drift=velblood+c1*(H0*dH);
		diffusion=c2;
%
		out1(ii1:ii2)=drift;
		out2(ii1:ii2)=diffusion;
		out3(ii1:ii2)=0;
		
	end
else
    
    switch(flag)
    case 'init'  
        out1 = t;
        
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
%::::::::::::::::::::::  DEFINE HERE THE SDE INITAL CONDITIONS  :::::::::::::::::::::::::::::::::
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%		out2=reshape((reshape(x0,3,1)*ones(1,NUMSIM)),1,NUMSIM*DIM);
		out2=x0;        
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: :::::::::::::::::::::::::::::::::
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        out3 = [];
        
        
    otherwise
        error(['Unknown flag ''' flag '''.']);
    end
end

