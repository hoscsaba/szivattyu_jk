function out = jk_rel_vel(z,C,S,geo)
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
  out.u= real(f);
  out.v=-imag(f);
end
