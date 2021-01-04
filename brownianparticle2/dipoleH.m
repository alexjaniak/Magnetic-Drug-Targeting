function H = dipoleH (m,r)
%DIPOLEH Calculate the magnetic field H of a diple in vacuum
%
%Input: 
%	m(1:3)	-- dipole moment
%	r(1:3)	-- relative location of the field point to the dipole
%Output:
%	H(1:3)	-- magnetic field at ponit r
%
%Record of revisions:
%	03/27/2009, Pengtao Yue, original code
%=====================

% Get the norm of r
rr=sqrt(sum(r(1:3).*r(1:3)));
% calculate H
H=(1/(4*pi))*(3*sum(m(1:3).*r(1:3))*r(1:3)/rr^5 - m(1:3)/rr^3);
