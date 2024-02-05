function out = jk_vel(z,C,S,geo,varargin)
  if length(z)>1
    error('!!!');
  end

  f=geo.Q_source/(2*pi)/z;
  for ll=1:geo.N_lapat
    for kk=1:geo.N_r-1
      z0=geo.x_c(kk,ll)+1i*geo.y_c(kk,ll);
      f=f-1i*C(kk)/(2*pi)/(z0-z)-S(kk)/(2*pi)/(z0-z);
    end
  end
  out.u= real(f)-abs(z)*geo.omega*sin(angle(z));
  out.v=-imag(f)+abs(z)*geo.omega*cos(angle(z));


  if nargin==5
    u_inf_x=varargin{1}(1);
    u_inf_y=varargin{1}(2);
    out.u=out.u+u_inf_x;
    out.v=out.v+u_inf_y;
  end
end
