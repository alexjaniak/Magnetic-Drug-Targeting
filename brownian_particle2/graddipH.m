function dh = graddipH(m,r)
%GRADDIPH Calculate the gradient of the H of a dipole
%
%Input: 
%	m(1:3)	-- diple moment
%	r(1:3)	-- relative location of the field point to the dipole
%Output:
%	dh(1:3,1:3)	-- grad H at r (dh(i,j)=\partial h(j)/ \partial x(j))
%Working:
%   mr      -- the inner product of m and r
%   rr      -- the norm of r
%   delta   -- the 3x3 identity matrix
%
%Record of revisions:
%	03/27/2009, Pengtao Yue, orignal code
%==================


dh=zeros(3,3);
mr=sum(m(1:3).*r(1:3));
rr=sqrt(sum(r(1:3).*r(1:3)));
delta=eye(3);

for ii=1:3
for jj=1:3
dh(ii,jj)=3*(m(ii)*r(jj)+mr*delta(ii,jj))/rr^5-15*mr*r(ii)*r(jj)/rr^7+ ...
          3*r(ii)*m(jj)/rr^5;
end
end
dh=dh/(4*pi);




