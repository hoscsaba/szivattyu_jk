function lapattervezes_v15
clear all, close all

%% Main inputs
Q_target=30/1000; H_target=40;
n=1440; g=9.81;
N_lapat=5; N_r=40; % min 4!!!
fname_prefix='jk_1';

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
geo.d_phi=pi*ones(1,N_r-1)/(N_r+1);
geo=jk_build_geo(geo);

%% Percentages of H/Q
A_C=0.11;
A_S=0*0.023;
[ff,geo]=obj(A_C,A_S,geo,0);

% fname=[fname_prefix,'_Q_',num2str(round(Q_target*3600)),'m3ph_H_',...
%     num2str(round(geo.H_target)),'m.mat']
% save(fname,"geo");

% geo.N_r-1
end
%% Defining the circulation
function C=get_C(A_C,xi,geo)
Int_C=0;
for i=1:length(xi)
    % xi=geo.loc_c(i)/geo.t_arclength(end);
    %C(i)=0.37;
    % Négyzetes
    % C(i)=A;
    % Parabolikus
    %C(i)=polyval(x,xi)*xi*(1-xi);
    %Elliptikus
    C(i)=0;
    if xi(i)<0.5
        C(i)=sin(acos(1-2*xi(i)));
    else
        C(i)=sin(acos(2*xi(i)-1));
    end
    % C
    if i==1
    Int_C=Int_C+C(i)*(xi(i));
    else
    Int_C=Int_C+C(i)*(xi(i)-xi(i-1));
    end
end
A=9.81*geo.H_target*2*pi*geo.b2/geo.omega/geo.N_lapat/geo.t_arclength(end)/Int_C;
%A=A_C;
C=C*A;
end

%% Defining the sources
function S=get_S(A_S,xi,geo)
Int_S=0;

for i=1:length(xi)
    S(i)=0;
    %xi=geo.loc_c(i)/geo.t_arclength(end);
    if xi(i)<0.3
        S(i)=sqrt(xi(i))*(0.33-xi(i))^2;
        if i==1
            Int_S=Int_S+S(i)*(xi(i));
        else
            Int_S=Int_S+S(i)*(xi(i)-xi(i-1));
        end
    elseif xi(i)>1-0.3
        S(i)=-sqrt(-xi(i)+1)*(0.33-1+xi(i))^2;
    else
        S(i)=0;
    end
end
A=A_S*geo.Q_target/Int_S/geo.N_lapat/geo.b2/geo.t_arclength(end);
S=S*A;
%S=S/4.1985e3*4.1985e3
%Int_S*geo.N_lapat*geo.b2*geo.t_arclength(end)
%geo.Q_target
% figure(124)
% plot(pxi,S)
end

%% Iteration
function [out,geo]=obj(A_C,A_S,geo,DO_PLOT)
xi=geo.loc_c/geo.t_arclength(end);

delta_d_phi=1e5;
iter=1;
ITER_MAX=50;
while (delta_d_phi>0.01) && (iter<ITER_MAX)
    C=get_C(A_C,xi,geo);
    S=get_S(A_S,xi,geo);
    d_phi_old=geo.d_phi;
    alpha=2*pi/geo.N_lapat/2;
    tmax=0.1;
    ode_options=odeset('Events',@(t,z) streamlineevent(t,z,C,S,geo),...
        'AbsTol',1e-5);
    [ts,xsys,te,xe,ie]=ode45(@(t,z) jk_streamlineode(t,z,C,S,geo),...
        [0 tmax],[geo.Db/2*cos(alpha);geo.Db/2*sin(alpha)],ode_options); %DE megoldó

%     figure(100)
%                    subplot(2,3,iter)
%     plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
%     title([num2str(iter),'. iteráció',])
% pause

    if length(ie)<geo.N_r
        fprintf('\n\n !!!!!!!!!!\n Last step: %g, tmax=%g\n',ts(end),tmax);
        fprintf('\n ==> increase tmax!');

        ie'

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

                figure(100)                
                plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')


                %pause
            else
                phi(ii)=atan2(xe(idx,2),xe(idx,1));
                if (phi(ii)<0)
                    phi(ii)=phi(ii)+2*pi;
                end
            end
        end

        aa=0.2;
        geo.d_phi=aa*diff(phi)+(1-aa)*d_phi_old;
        delta_d_phi=norm(geo.d_phi-d_phi_old);
        geo=jk_build_geo(geo);

        %DO_PLOT=1;

        [QQ,HH,veldata,geo]=jk_kompl_pot(C,S,geo,DO_PLOT);
        geo.QQ=QQ;
        geo.HH=HH;
        geo.veldata=veldata;

        err1=(geo.H_target-HH)^2;
        err2=0*std(geo.veldata.c_k_u_vec);
    end
    fprintf('\n\t iter #%d, norm(d_phi change) = %5.3e',iter,delta_d_phi);

    iter=iter+1;

    figure(100)

    plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
    title(['Q=',num2str(round(QQ*1000)),...
        'l/s, H=',num2str(round(10*HH)/10),'m, iter #',num2str(iter)])
%pause
if iter==ITER_MAX
    warning("Blade iteration not converging!");
    pause
end
end

out=0*err1+5*err2;
geo.C=get_C(A_C,xi,geo);
geo.S=get_S(A_S,xi,geo);
        figure(100)
        subplot(2,3,[1,2,4,5])
        plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
        subplot(2,3,3)
        xx=linspace(0,1);
        plot(xx,get_C(A_C,xx,geo))
        subplot(2,3,6)
        plot(xx,get_S(A_S,xx,geo))

%% Send the geometry to postprocessing
geo=jk_postprocess(geo);

%% Save to CFX
jk_postprocess(geo);
%save_to_CFX(geo);

end




function [val,ter,dir]=streamlineevent(t,z,C,S,geo)

Ract=sqrt(z(1)^2+z(2)^2);
R_e=linspace(geo.Db,geo.D2,geo.N_r)/2;
val=Ract-R_e;
ter=zeros(size(val)); ter(end)=1;
dir=zeros(size(val));
end
