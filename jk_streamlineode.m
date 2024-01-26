function dydx=jk_streamlineode(t,z,C,S,geo)
tmp=jk_vel(z(1)+1i*z(2),C,S,geo);
dydx=[tmp.u;tmp.v];
end