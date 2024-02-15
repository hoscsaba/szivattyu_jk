function geo=jk_postprocess(geo)
    % figure(123),clf, hold on
    % subplot(2,3,[1,2,4,5])
    % 
    % plot(geo.x_g,geo.y_g,'k','LineWidth',1)
    % hold on
    % r_b=geo.Db/2;
    % r_k=geo.D2/2;
    % 
    % fi_ini=pi*linspace(0,360,100)/180;
    % ode_options=odeset('Events',@(t,z) jk_streamlineevent(t,z,geo.C,geo.S,geo));
    % for ii=1:length(fi_ini)
    %     xini=0.8*r_b*cos(fi_ini(ii));
    %     yini=0.8*r_b*sin(fi_ini(ii));
    %     [ts,xsys]=ode45(@(t,z) jk_streamlineode(t,z,geo.C,geo.S,geo),...
    %         [0 0.05],[xini;yini],ode_options); %DE megoldó
    % 
    %     plot(xsys(:,1),xsys(:,2),'b');
    % end
    % 
    % hold off
    % 
    % title(['N_r=',num2str(geo.N_r),', Q=',num2str(round(geo.QQ*3600)),...
    %     ' m3/h, H=',num2str(round(geo.HH*100)/100),'m'])
    % axylim=r_k*1.5;
    % xlim([-axylim,axylim]), ylim([-axylim,axylim]);
    % axis("equal");
    % 
    % drawnow


% function jk_lapat_deformacio
% clear all, close all
% 
% global L a EI geo w_lapat
% 
% fname='jk_1p0_Q_170m3ph_H_87m.mat';
% rho_lapat=7860; h_lapat=3e-3;
% geo=jk_cp(fname,1,rho_lapat,h_lapat);
% 
% L=geo.t_arclength(end); a=L/3; EI=2.8; w_lapat=geo.b2;
% 
% x=linspace(0,1);
% solinit=bvpinit(x,@guess);
% sol=bvp5c(@tarto_ode,@tarto_bv,solinit);
% 
% figure(1)
% subplot(2,2,1)
% plot(L*geo.x_dp,geo.dp,'r-',...
%      L*geo.x_dp,geo.pcp,'b-')
%      xlabel('Lapát ívhossz'), ylabel('p_{terhelés}'), legend('\Deltap','p_{cp}')
%      xlim([0,0.534])
% 
% 
% subplot(2,2,3)
% plot(a*sol.x,sol.y(1,:),...
%     a+(L-2*a)*sol.x,sol.y(5,:),...
%     (L-a)+a*sol.x,sol.y(9,:),...
%     a,0,'ro',L-a,0,'ro')
%     xlabel('Lapát ívhossz'), ylabel('deformáció')
%     xlim([0,0.534])
% 
% 
% %% deformalt lapatalak
% def_x=[a*sol.x(1:end-1), a+(L-2*a)*sol.x(1:end-1), (L-a)+a*sol.x];
% def_x=def_x/def_x(end);
% def_y=[sol.y(1,1:end-1), sol.y(5,1:end-1),sol.y(9,:)];
%  for i=1:length(geo.x_c(1,:))
%      for np=1:length(geo.x_c(:,i))
%         n_v=[geo.n_x(np); geo.n_y(np)];
%         ds=geo.t_arclength(np)/geo.t_arclength(end);
%         def=interp1(def_x,def_y,ds);
%         x_d(np,i)=geo.x_c(np,i)+def*n_v(1);
%         y_d(np,i)=geo.y_c(np,i)+def*n_v(2);         
%      end
%  end
% 
% subplot(2,2,[2,4])
% figure(1)
%  for i=1:1%length(x_d(1,:))
%     plot(x_d(1:end-1,i),y_d(1:end-1,i),'b',LineWidth=1), hold on      
%  end
%  for i=1:length(geo.x_g(1,:))
%     plot(geo.x_g(:,i),geo.y_g(:,i),'k'), hold on      
%  end
% axis equal
% legend('deformált lapátalak','eredeti lapátalak')
% xlim([-0.3,0.3])
% 
% figure(8)
%  for i=1:1
%     plot(x_d(1:end-1,i),y_d(1:end-1,i),'b',LineWidth=1), hold on      
%  end
%  for i=1:1
%     plot(geo.x_g(:,i),geo.y_g(:,i),'k'), hold on      
%  end
% axis equal
% legend('deformált lapátalak','eredeti lapátalak')
% xlim([-0.3,0.3])
% 
% figure(10)
%  for i=1:1%length(x_d(1,:))
%     plot(x_d(1:end-1,i),y_d(1:end-1,i),'b',LineWidth=1), hold on      
%  end
%  for i=1:length(geo.x_g(1,:))
%     plot(geo.x_g(:,i),geo.y_g(:,i),'k'), hold on      
%  end
% axis equal
% %legend('deformált lapátalak','eredeti lapátalak')
% xlim([-0.3,0.3])
%      pbaspect([1 1 1])
% 
% end
% 
% %% Inicializálás (nem lényeges, csak kell vmi)
% function y=guess(x)
% for i=1:12, y(i,1)=0; end
% end
% 
% %% Terhelés, ezt kell majd megadnod!
% function out = q(x)
% global L a geo
% 
% xx=x/L;
% out_p =-interp1(geo.x_dp,geo.dp,xx);
% out_cp=-interp1(geo.x_dp,geo.pcp,xx);
% 
% % Aszimmetrikus négyzetes terhelés
% %out=-x^2;
% 
% % Szimmetrikus négyzetes terhelés
% %out=-(x-L/2)^2;
% 
% % Terhelés csak középen
% %if x<a
% %    out=0;
% %elseif x<L-a
% %    out=-1;
% %else
% %    out=0;
% %end
% 
% out=(out_p+out_cp)*geo.b2*30;
% %out=(out_p+out_cp)*geo.b2;
% end
% 
% %% ODE
% function dzdx = tarto_ode(x,y)
% global L a EI
% % 1.szakasz       2.szakasz    3.szakasz
% %  =============x=============x===========
% % |             |             |           |
% % |<----------->|<----------->|<--------->|
% %        a                          a
% % |<------------------------------------->|
% %                      L
% % z(1): 1. szakasz y
% % z(2): 1. szakasz dy/dx
% % z(3): 1. szakasz d2y/dx2
% % z(4): 1. szakasz y3y/dx3
% 
% % z(5): 2. szakasz y
% % z(6): 2. szakasz dy/dx
% % z(7): 2. szakasz d2y/dx2
% % z(8): 2. szakasz y3y/dx3
% 
% % z(9):  3. szakasz y
% % z(10): 3. szakasz dy/dx
% % z(11): 3. szakasz d2y/dx2
% % z(12): 3. szakasz y3y/dx3
% 
% % x=0...1
% 
% % 1. szakasz
% xx1=x*a;
% dzdx(1) = a*y(2); % dz1/dx=z2
% dzdx(2) = a*y(3); 
% dzdx(3) = a*y(4); 
% dzdx(4) = a*q(xx1)/EI; 
% 
% % 2. szakasz
% xx2=a+(L-2*a)*x;
% dzdx(5) = (L-2*a)*y(6);
% dzdx(6) = (L-2*a)*y(7); 
% dzdx(7) = (L-2*a)*y(8); 
% dzdx(8) = (L-2*a)*q(xx2)/EI; 
% 
% % 3. szakasz
% xx3=(L-a)+x*a;
% dzdx(9)  = a*y(10);
% dzdx(10) = a*y(11); 
% dzdx(11) = a*y(12); 
% dzdx(12) = a*q(xx3)/EI; 
% end
% 
% %% Peremfeltételek
% function out=tarto_bv(ya,yb)
% % bal vége szabad, terheletlen y''=0 és y'''=0 @ 0
% out(1)=ya(3);
% out(2)=ya(4);
% % jobb vége szabad, terheletlen  y''=0 és y'''=0 @ 1
% out(3)=yb(11);
% out(4)=yb(12);
% % bal megtámasztás
% out(5)=yb(1);
% % jobb megtámasztás
% out(6)=ya(9);
% % bal csatlakozás
% out(7)=yb(1)-ya(5);
% out(8)=yb(2)-ya(6);
% out(9)=yb(3)-ya(7);
% % jobb csatlakozás
% out(10)=yb(5)-ya(9);
% out(11)=yb(6)-ya(10);
% out(12)=yb(7)-ya(11);
% end
% 
% 
% 
% function geo=jk_cp(fname,DO_PLOT,rho_lapat,h_lapat)
% 
% %fname='jk_0p8_Q_136m3ph_H_104m.mat';
% load(fname);
% %load('jk_0p8_Q_136m3ph_H_104m.mat');
% C=geo.C;
% geo.rho=1000;
% geo.rho_lapat=rho_lapat;
% geo.h_lapat=h_lapat;
% 
% tmax=0.1;
% ode_options=odeset('Events',@(t,z) streamlineevent(t,z,C,geo),'AbsTol',1e-5);
% R0=0.9*geo.Db/2;
% alpha=1*pi/180;
% [t_s,x_s,te,xe,ie]=ode45(@(t,z) streamlineode(t,z,C,geo),...
%     [0 tmax],[R0*cos(alpha);R0*sin(alpha)],ode_options); %DE megoldó
% 
% alpha=-40*pi/180;
% [t_n,x_n,te,xe,ie]=ode45(@(t,z) streamlineode(t,z,C,geo),...
%     [0 tmax],[R0*cos(alpha);R0*sin(alpha)],ode_options); %DE megoldó
% 
% 
% if DO_PLOT==1
    figure(5)
    subplot(2,2,[1,3])
    plot(geo.x_g,geo.y_g,'k','LineWidth',1)
    hold on
    r_b=geo.Db/2;
    r_k=geo.D2/2;

    ode_options=odeset('Events',@(t,z) jk_streamlineevent2(t,z,geo.C,geo.S,geo));
    alpha=-27*pi/180;
    [t_s,xsys_s,te,xe,ie]=ode45(@(t,z) jk_streamlineode(t,z,geo.C,geo.S,geo),...
        [0 0.05],[0.5*r_b*cos(alpha);0.5*r_b*sin(alpha)],ode_options); %DE megoldó
    plot(xsys_s(:,1),xsys_s(:,2),'r','LineWidth',2);

    alpha=-39*pi/180;
    [t_n,xsys_n,te,xe,ie]=ode45(@(t,z) jk_streamlineode(t,z,geo.C,geo.S,geo),...
        [0 0.05],[0.5*r_b*cos(alpha);0.5*r_b*sin(alpha)],ode_options); %DE megoldó
    plot(xsys_n(:,1),xsys_n(:,2),'g','LineWidth',2);

    [s_s,r_s,w_s,p_s,p_cp]=jk_get_w(xsys_s,geo);
    [s_n,r_n,w_n,p_n,p_cp]=jk_get_w(xsys_n,geo);

    fi_ini=pi*linspace(0,360,100)/180;
    ode_options=odeset('Events',@(t,z) jk_streamlineevent2(t,z,geo.C,geo.S,geo));
    for ii=1:length(fi_ini)
        xini=0.8*r_b*cos(fi_ini(ii));
        yini=0.8*r_b*sin(fi_ini(ii));
        [ts,xsys]=ode45(@(t,z) jk_streamlineode(t,z,geo.C,geo.S,geo),...
            [0 0.05],[xini;yini],ode_options); %DE megoldó

        plot(xsys(:,1),xsys(:,2),'b');
    end



    %% Interpolacio
    l=linspace(0,geo.t_arclength(end),101);
    w_s_ip=interp1(s_s,w_s,l);
    w_n_ip=interp1(s_n,w_n,l);
    delta_w_ip=w_s_ip-w_n_ip;

    %Gamma_lapat_testing=trapz(l,delta_w_ip)


    %Gamma_lapat_testing=Gamma_lapat_testing*geo.t_arclength(end)/100

    hold off

    subplot(2,2,2)
    u2=geo.D2/2*geo.omega;
    plot(2*r_s/geo.Db,w_s/u2,'r'), hold on
    plot(2*r_n/geo.Db,w_n/u2,'b')
    xlim([1,max(2*r_s/geo.Db)])
    xlabel('r/D_1'), ylabel('w/u_2'), legend('szívott oldal','nyomott oldal')

    subplot(2,2,4)
    plot(s_s,p_s,'r'), hold on
    plot(s_n,p_n,'b')
    xlabel('Áramvonal ívhossz'), ylabel('p'), legend('szívott oldal','nyomott oldal')
    
% 
% %% Interpolacio
% xx1=linspace(0,1);
% pps=interp1(s_s/s_s(end),p_s,xx1);
% ppn=interp1(s_n/s_n(end),p_n,xx1);
% pcp=interp1(s_n/s_n(end),p_cp,xx1);
% 
% geo.x_dp=xx1;
% geo.dp=ppn-pps;
% geo.pcp=pcp;
% 

    
end