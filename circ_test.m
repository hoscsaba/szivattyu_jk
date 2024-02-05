function circ_test
  clear all, close all, clc

  global geo Ux Uy

  Ux=1; Uy=0;


  alpha=5; L=1;
  geo.N_lapat=1;
  geo.N_r=10;
  for j=1:5
    geo.x_c=linspace(0,L,geo.N_r)';
    geo.y_c=linspace(0,-L*tand(alpha),geo.N_r)';
    geo.D2=10.;

    geo.Q_source=0;
    geo.omega=0;

    for i=1:geo.N_r-1
      geo.xp(i)=(geo.x_c(i)+geo.x_c(i+1))/2;
      geo.yp(i)=(geo.y_c(i)+geo.y_c(i+1))/2;
      geo.nxp(i)=(geo.y_c(i)-geo.y_c(i+1))/2;
      geo.nyp(i)=-(geo.x_c(i)-geo.x_c(i+1))/2;
      tmp=sqrt(geo.nxp(i)^2+geo.nyp(i)^2);
      geo.nxp(i)=geo.nxp(i)/tmp;
      geo.nyp(i)=geo.nyp(i)/tmp;
    end

    Cini=0.1*ones(1,geo.N_r);
    SS=0.0;
    for i=1:geo.N_r
      if i<geo.N_r/3
        geo.S(i)=SS;
      elseif i<2*geo.N_r/3
        geo.S(i)=0;
      else
        geo.S(i)=-SS;
      end
    end

    C=fsolve(@BC,Cini(1:end-1))
    geo.C=[C 0];

    Np=20;
    yini_v=linspace(-0.2,0.1,Np);
    ode_options=odeset('Events',@(t,z) jk_streamlineevent(t,z,geo.C,geo.S,geo));

    figure(j)
    subplot(2,1,1)
    for i=1:Np
      [ts,xsys]=ode45(@(t,z) jk_streamlineode(t,z,geo.C,geo.S,geo,[Ux Uy]),...
        [0 2],[-0.4;yini_v(i)],ode_options); %DE megoldó
      plot(xsys(:,1),xsys(:,2),'b'), hold on
    end
    plot(geo.x_c,geo.y_c,'r*'), hold on
    plot(geo.xp,geo.yp,'bo'), hold on
    mul=0.1; mul2=0.1;
    for i=1:geo.N_r-1
      plot([geo.xp(i) geo.xp(i)+mul*geo.nxp(i)],[geo.yp(i) geo.yp(i)+mul*geo.nyp(i)],'k')
      tmp=jk_vel(geo.xp(i)+geo.yp(i)*1i,geo.C,geo.S,geo,[Ux Uy]);
      plot([geo.xp(i) geo.xp(i)+mul2*tmp.u],[geo.yp(i) geo.yp(i)+mul2*tmp.v],'r')
    end
    hold off
    %axis equal
    title(["N_r=",num2str(geo.N_r),", sum(C)=",num2str(sum(geo.C))]);
    subplot(2,1,2)
    plot(C,'-o');
drawnow
    geo.N_r=geo.N_r*2;
  end
end

function out = BC(x)
  global geo Ux Uy

  for i=1:geo.N_r-1
    z   = geo.xp(i)+1i*geo.yp(i);
    tmp = jk_vel(z,[x 0],geo.S,geo,[Ux,Uy]);
    v   = [tmp.u tmp.v];
    n   = [geo.nxp(i) geo.nyp(i)];
    out(i) = dot(n,v);
  end
end
