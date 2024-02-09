function lapattervezes_v13
  clear all, close all

  %% Main inputs
  Q_target=20/60; H_target=50;
  n=1450; g=9.81;
  N_lapat=5; N_r=40; % min 4!!!

  %% Calculating the main parameters of the pump
  nq=n*Q_target^0.5/H_target^0.75; psi=(300/(300+nq))^(9/4);
  u2=sqrt(2*g*H_target/psi); D2=u2/(pi*n/60); Db=D2*0.4;
  omega=2*pi*n/60;
  u1=Db*pi*n/60; u2=D2*pi*n/60;
  geo.u2=u2;
  epszilon=0.0188*nq^(2/3); c1=epszilon*sqrt(2*g*H_target);
  beta1=atan(c1/u1);
  c2u=H_target*g/u2; c2m=0.1011*sqrt(2*g*H_target);
  geo.beta2=atan(c2m/(u2-c2u));
  b2=Q_target/(D2*pi*c2m);
  b1=Q_target/(Db*pi*c1);

  fname_prefix='jk_A';

  fprintf("\n Q  = %5.1f l/s",Q_target*60);
  fprintf("\n H  = %5.1f m",H_target);
  fprintf("\n nq = %5.3f",nq);
  fprintf("\n n  = %5.3f rpm",n);
  fprintf("\n D1 = %5.3f mm",Db*1000);
  fprintf("\n D2 = %5.3f mm",D2*1000);
  fprintf("\n b2 = %5.3f mm",b2*1000);
  pause
  %Declaring the geo structure
  geo.Q_target=Q_target;
  geo.Q_source=Q_target/b2;
  geo.N_lapat=N_lapat;
  geo.N_r=N_r;
  geo.b2=b2;
  geo.Db=Db;
  geo.D2=D2;
  geo.beta1=beta1;
  geo.omega=omega;
  geo.H_target=H_target;
  geo.Gamma_lapat_elm=9.81*H_target*2*pi/N_lapat/omega;

  %% Initial geometry, for temporary purposes
  geo.d_phi=2*pi*ones(1,N_r-1)/(N_r+1);
  geo=jk_build_geo(geo);

  %% Run a series of computations
  [ff,geo]=obj(0.15,geo,0);

  save_to_CFX(geo);

end
%% Calculating the circulation
function C=get_C(A,geo)
  for i=1:geo.N_r-1
    xi=geo.loc_c(i)/geo.t_arclength(end);
    if xi<0.5
      C(i)=A*sin(acos(1-2*xi));
    else
      C(i)=A*sin(acos(2*xi-1));
    end
  end
end

%% Calculating the sources
function S=get_S(A,geo)
  A=0;
  for i=1:geo.N_r-1
    xi=geo.loc_c(i)/geo.t_arclength(end);
    if xi<0.3
      S(i)=-4*sqrt(xi)*(0.33-xi)^2;
    elseif xi>1-0.3
      S(i)=4*sqrt(-xi+1)*(0.33-1+xi)^2;
    else
      S(i)=0;
    end
    pxi(i)=xi;
  end
  S=A*S;
end

%% Error function
function [out,geo]=obj(A,geo,DO_PLOT)

  C=get_C(A,geo);
  S=get_S(0,geo);

  delta_d_phi=1e5;
  iter=1;
  ITER_MAX=20;
  while (delta_d_phi>0.01) && (iter<ITER_MAX)

    d_phi_old=geo.d_phi;
    alpha=2*pi/geo.N_lapat/2;
    tmax=0.2;
    ode_options=odeset('Events',@(t,z) streamlineevent(t,z,C,S,geo),...
      'AbsTol',1e-5);
    [ts,xsys,te,xe,ie]=ode45(@(t,z) streamlineode(t,z,C,S,geo),...
      [0 tmax],0.95*[geo.Db/2*cos(alpha);geo.Db/2*sin(alpha)],ode_options); %DE megoldó

%    figure(100)
 %   subplot(2,3,iter)
  %  plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
  %  title([num2str(iter),'. iteráció',])
  %  pause

    if length(ie)<geo.N_r
      fprintf('\n\n !!!!!!!!!!\n Last step: %g, tmax=%g\n',ts(end),tmax);
      fprintf('\n ==> increase tmax!');

      ie'
      figure(100)
      plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r',xe(:,1),xe(:,2),'ro')

      %pause
      err1=1e5;
      err2=1e5;
      HH=0;
    else
      for ii=1:geo.N_r
        idx=find(ie==ii,1,'first');
        if isempty(idx)
          phi(ii)=15*pi/180;
          warning('!!!!!!!!!!\n dr not found, replacing by 15 degrees');
          ie'
          pause
          figure(100)                
          plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r',xe(:,1),ye(:,2),'ro')


          %pause
        else
          phi(ii)=atan2(xe(idx,2),xe(idx,1));
          %if (phi(ii)<0)
          %  phi(ii)=phi(ii)+2*pi;
          %end
        end
      end
      aa=1;
      geo.d_phi=aa*diff(phi)+(1-aa)*d_phi_old;
      delta_d_phi=norm(geo.d_phi-d_phi_old);
      geo=jk_build_geo(geo);

      DO_PLOT=0;

      [QQ,HH,veldata,geo]=jk_kompl_pot(C,S,geo,DO_PLOT);

      err1=(geo.H_target-HH)^2;
      err2=0*std(veldata.c_k_u_vec);
    end
    fprintf('\n\t iter #%d, norm(d_phi change) = %5.3e',iter,delta_d_phi);

    iter=iter+1;

    figure(100)

    plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
    title(['Q=',num2str(round(QQ*60)),...
      'l/s, H=',num2str(round(10*HH)/10),'m, iter #',num2str(iter)])
  %pause
  if iter==ITER_MAX
    warning("Blade iteration not converging!");
    pause
  end
  end
  % fprintf('\n Q=%5.1f m3/h, H=%5.1f / %5.1f m, err=[ %5.3f, %5.3f], x=[',...
  %     QQ*3600,geo.H_target,HH,err1,err2);
  % for ii=1:length(x)
  %     fprintf(' %5.3f ',x(ii));
  % end
  % fprintf("]")

  out=err1+err2;
end

function dydx=streamlineode(t,z,C,S,geo)
  tmp=jk_vel(z(1)+1i*z(2),C,S,geo);
  dydx=[tmp.u;tmp.v];
end

function [val,ter,dir]=streamlineevent(t,z,C,S,geo)

  Ract=sqrt(z(1)^2+z(2)^2);
  R_e=linspace(geo.Db,geo.D2,geo.N_r)/2;
  val=Ract-R_e;
  ter=zeros(size(val)); ter(end)=1;
  dir=zeros(size(val));
end
