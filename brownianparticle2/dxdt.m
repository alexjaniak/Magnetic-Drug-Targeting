function dx=dxdt(tt,xx)
%DXDT Returns velocity to be sent to odesolvers
%
%   xx'(tt)=dxdt
%
%Input:
%	tt	-- time, independent variable
%	xx(3,1)	-- particle position, column vector
%
%Output:
%	dx(3,1) -- xx'(tt)
%
%Global:
%	m_dpl(3) -- magnetization of the magnet
%	mu0	-- permeability of vacuum
%	xdipole -- position of the magnet
%	rad_particle -- radius of particle
%	chi	-- susceptibility of the ferro particle
%	rad_pipe -- radius of the pipe centered along x-axis
%	viscosity -- viscosity of the flowing liquid
%	vmax	-- mean max flow velocity in the pipe
%	vosc	-- amplitude of velocity oscillation
%	omega	-- angular frequency of oscillation
%Record of revisions:
%	05/27/2009, Pengtao Yue, original code
%	06/12/2009, Pengtao Yue, add oscillatory background flow
%==========================
global m_dpl chi mu0 viscosity rad_pipe rad_particle xdipole vmax vosc omega

x=reshape(xx,1,3);
vbkg=zeros(1,3);
vbkg(1)=poiseuille(x(1:3),rad_pipe,vmax,vosc,omega,tt);

vslp=slipvel(m_dpl,mu0,x-xdipole,chi,rad_particle,viscosity);

dx=reshape(vbkg+vslp,3,1);


