function dydx=jk_streamlineode2(t,z,C,S,geo,varargin)
if nargin==6
  tmp=jk_vel2(z(1)+1i*z(2),C,S,geo,varargin{1});
else
  tmp=jk_vel2(z(1)+1i*z(2),C,S,geo);
end
dydx=[tmp.u;tmp.v];
end
