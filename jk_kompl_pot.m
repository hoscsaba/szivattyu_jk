function [QQ,HH,veldata,geo]=jk_kompl_pot(C,geo,DO_PLOT)

[QQ,HH,veldata.phi_ker_vec,veldata.c_k_u_vec,veldata.c_k_m_vec]=Get_QH(C,geo);

veldata.Gamma_lapat=Get_Circulation(geo,C);

if DO_PLOT==1
    figure(123),clf, hold on
    subplot(2,3,[1,2,4,5])

    plot(geo.x_g,geo.y_g,'k','LineWidth',1)
    hold on
    RR=1.01;
    mul=1/200;
    r_b=geo.Db/2;
    r_k=geo.D2/2;
    phi_vec=linspace(0,2*pi,7);
    for ii=1:length(phi_vec)
        phi_ker=phi_vec(ii);
        n_k=[cos(phi_ker);sin(phi_ker)];
        t_k=[sin(phi_ker);-cos(phi_ker)];
        x_r_k=RR*r_k*cos(phi_ker);
        y_r_k=RR*r_k*sin(phi_ker);
        velo=vel(x_r_k+y_r_k*1i,C,geo);
        w_k=[velo.u;velo.v];
        u_k=RR*r_k*geo.omega*t_k;
        %c_k=w_k+u_k;
        dv1=mul*u_k;
        dv2=mul*w_k;
        dv3=mul*(u_k+w_k);

        h1=annotation('arrow');
        % [x0,y1,dx,fy]
        set(h1,'parent', gca, 'position', [x_r_k, y_r_k,  dv1(1), dv1(2)],...
            'Color',[0 0 1], 'LineWidth', 1);

        h2=annotation('arrow');
        set(h2,'parent', gca, 'position', [x_r_k+dv1(1), y_r_k+dv1(2),  dv2(1), dv2(2)],...
            'Color',[0 1 0], 'LineWidth', 1);

        h3=annotation('arrow');
        set(h3,'parent', gca, 'position', [x_r_k, y_r_k,  dv3(1), dv3(2)],...
            'Color',[1 0 0], 'LineWidth', 1);

    end

    fi_ini=linspace(0,2*pi/geo.N_lapat,30);
    ode_options=odeset('Events',@(t,z) streamlineevent(t,z,C,geo));
    for ii=1:length(fi_ini)
        xini=0.99*r_b*cos(fi_ini(ii));
        yini=0.99*r_b*sin(fi_ini(ii));
        [ts,xsys]=ode23(@(t,z) streamlineode(t,z,C,geo),...
            [0 0.1],[xini;yini],ode_options); %DE megoldó

        plot(xsys(:,1),xsys(:,2),'r');
    end

    plot(geo.x_g,geo.y_g,'ro')
    plot(geo.x_c,geo.y_c,'k^')
    plot(geo.x_v,geo.y_v,'b*')
    mul2=0.2;
    for ii=1:geo.N_r-1
        h1=annotation('arrow');
        set(h1,'parent', gca, 'position',...
            [geo.x_v(ii), geo.y_v(ii),  mul2*geo.n_x(ii), mul2*geo.n_y(ii)],...
            'Color',[0 0 1], 'LineWidth', 1);
    end

    mul3=0.01;
    for ii=1:geo.N_r-1
        dx=1e-3/10;
        z=geo.x_v(ii)+1i*geo.y_v(ii);
        vf=vel(z + (geo.n_x(ii)+1i*geo.n_y(ii))*dx,C,geo);
        va=vel(z - (geo.n_x(ii)+1i*geo.n_y(ii))*dx,C,geo);
        h1=annotation('arrow');
        set(h1,'parent', gca, 'position',...
            [geo.x_v(ii)+dx*geo.n_x(ii), geo.y_v(ii)+dx*geo.n_y(ii),...
            mul3*vf.u, mul3*vf.v],...
            'Color',[1 0 0], 'LineWidth', 1);
        h2=annotation('arrow');
        set(h2,'parent', gca, 'position',...
            [geo.x_v(ii)-dx*geo.n_x(ii), geo.y_v(ii)-dx*geo.n_y(ii),...
            mul3*va.u, mul3*va.v],...
            'Color',[0 1 0], 'LineWidth', 1);
    end

    hold off

    title(['N_r=',num2str(geo.N_r),', Q=',num2str(round(QQ*3600)),...
        ' m3/h, H=',num2str(round(HH*100)/100),'m'])
    xylim=r_k*2;
    xlim([-xylim,xylim]), ylim([-xylim,xylim]);
    axis("equal");

    subplot(2,3,3)
    plot(veldata.phi_ker_vec,veldata.c_k_u_vec,...
        veldata.phi_ker_vec,veldata.c_k_m_vec)
    xlabel('\phi'), ylabel('c'), legend('c_u','c_m')

    subplot(2,3,6)
    %xx=linspace(0,1);
    plot(geo.loc_c,C,'o-')
    %plot(xx*t_arclength(end),Get_C_at_s(coeffs_opt,xx),'--');
    xlabel('lapát ívhossz'), ylabel('C')
    drawnow
    
end

end

function out=Get_Circulation(geo,C)

dx=1e-2;

ds=geo.t_arclength(end)/(geo.N_r-1);
%fprintf('\n');
for ii=1:geo.N_r-1
    z=geo.x_v(ii)+1i*geo.y_v(ii);
    vf=vel(z + (geo.n_x(ii)+1i*geo.n_y(ii))*dx,C,geo);
    va=vel(z - (geo.n_x(ii)+1i*geo.n_y(ii))*dx,C,geo);
    %fprintf('\n ii=%g, ds=%5.3f, vf=%5.3e, va=%5.3e, dv=%5.3e',...
    %   ii,ds,sqrt(vf.u^2+vf.v^2),sqrt(va.u^2+va.v^2),...
    %    sqrt(vf.u^2+vf.v^2)-sqrt(va.u^2+va.v^2));
    gg(ii)=ds*(sqrt(vf.u^2+vf.v^2)-sqrt(va.u^2+va.v^2));    
end
out=ds*sum(gg);

end

function dydx=streamlineode(t,z,C,geo)
tmp=vel(z(1)+1i*z(2),C,geo);
dydx=[tmp.u;tmp.v];
end

function [val,ter,dir]=streamlineevent(t,z,C,geo)
val=z(1)^2+z(2)^2-(geo.D2/2)^2;
ter=1;
dir=0;
end

function [Q,H,phi_ker_vec,c_k_u_vec,c_k_m_vec]=Get_QH(C,geo)

RR=1.01;
r_k=geo.D2/2;
phi_ker_vec=linspace(0,2*pi,50*geo.N_lapat);
for i_phi=1:length(phi_ker_vec)
    phi_ker=phi_ker_vec(i_phi);
    n_k=[cos(phi_ker);sin(phi_ker)];
    t_k=[sin(phi_ker);-cos(phi_ker)];
    x_r_k=RR*r_k*cos(phi_ker);
    y_r_k=RR*r_k*sin(phi_ker);
    velo=vel(x_r_k+y_r_k*1i,C, geo);
    w_k=[velo.u;velo.v];
    u_k=RR*r_k*geo.omega*t_k;
    c_k=w_k+u_k;
    c_k_u_vec(i_phi)=dot(c_k,t_k);
    c_k_m_vec(i_phi)=dot(c_k,n_k);
end

Q=mean(c_k_m_vec)*2*r_k*pi*geo.b2;
H=mean(c_k_u_vec)*norm(u_k)/9.81;
end
