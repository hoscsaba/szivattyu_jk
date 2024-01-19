function lapattervezes_v9

clear all, close all

Q_target=170/3600; H_target=87;
n=1450; g=9.81;
N_lapat=5; N_r=20; % min 4!!!

nq=n*Q_target^0.5/H_target^0.75; psi=(300/(300+nq))^(9/4);
u2=sqrt(2*g*H_target/psi); D2=u2/(pi*n/60); Db=D2*0.4;
omega=2*pi*n/60;

u1=Db*pi*n/60; u2=D2*pi*n/60;
epszilon=0.0188*nq^(2/3); c1=epszilon*sqrt(2*g*H_target);
beta1=atan(c1/u1); c2u=H_target*g/u2; c2m=0.1011*sqrt(2*g*H_target);
b2=Q_target/(D2*pi*c2m); b1=Q_target/(Db*pi*c1);

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


%% Initial geometry, fro temporary purposes
geo.d_phi=3/2*pi*ones(1,N_r-1)/(N_r+1);
geo=jk_build_geo2(geo);

%% Unloaded blade
C=zeros(N_r-1,1);
ode_options=odeset('Events',@(t,z) streamlineevent(t,z,C,geo));
[ts,xsys,te,xe,ie]=ode23(@(t,z) streamlineode(t,z,C,geo),...
    [0 1],[geo.Db/2;0],ode_options); %DE megoldó
%plot(xsys(:,1),xsys(:,2),'r',xe(:,1),xe(:,2),'bo');
for ii=1:length(xe(:,1))
    phi(ii)=atan2(xe(ii,2),xe(ii,1));
    if (phi(ii)<0)
        phi(ii)=phi(ii)+2*pi;
    end
end
geo.d_phi=diff(phi);

geo=jk_build_geo2(geo);
[QQ,HH,veldata,geo]=jk_kompl_pot(C*0,geo,1);
pause
C=solve_for_C(geo,0*ones(N_r-1,1),0);
[QQ,HH,veldata,geo]=jk_kompl_pot(C,geo,1);
pause

d_fi_min=0*pi/180;
d_fi_max=30*pi/180;
lb=ones(N_r-1,1)*d_fi_min;
ub=ones(N_r-1,1)*d_fi_max;

%options = optimoptions(@fmincon,'OptimalityTolerance',0.1);
%d_phi_opt=fmincon(@(x)obj(x,geo,0),geo.d_phi,[],[],[],[],lb,ub,[],options);

options = optimoptions('ga','PlotFcn', @gaplotbestf);
d_phi_opt = ga(@(x)obj(x,geo,0),length(geo.d_phi),[],[],[],[],lb,ub,[],[],options);

for i=1:length(d_phi_opt)
    fprintf("\n d_phi(%2d)= %5.3f -> %5.3f fok",...
        i,geo.d_phi(i)*180/pi,d_phi_opt(i)*180/pi);
end

%geo.d_phi=[geo.d_phi1, d_phi_opt];
geo=jk_build_geo2(geo);
C=solve_for_C(geo,ones(N_r-1,1),0);
jk_kompl_pot(C,geo,1);
end

function out=obj(x,geo,DO_PLOT)

geo.d_phi=x;
geo=jk_build_geo2(geo);
C=solve_for_C(geo,ones(geo.N_r-1,1),DO_PLOT);
[Q_act,H_act,veldata,~]=jk_kompl_pot(C,geo,DO_PLOT);

fprintf('\n Q=%5.1f m3/h, H=%5.1f (target), %5.1f (actual) m, C=%5.3f...%5.3f',...
    Q_act*3600,geo.H_target,H_act,min(C),max(C));

penalty=0.;
if min(C)*max(C)<0
    penalty=10000;
end

out=(geo.H_target-H_act)^2+std(veldata.c_k_u_vec)+penalty;
end

function out = Get_df1(geo)

Rb=geo.Db/2;
dr=(geo.D2-geo.Db)/geo.N_r/2;

x=dr*2;
dfi1=-60*pi/180;
for i=1:10
    z=(Rb^2+(Rb+dr)^2-2*Rb*(Rb+dr)*cos(dfi1)).^0.5;
    y=geo.Db+dr-sqrt(x^2+Rb^2);
    x=sqrt(y^2+2*x*z*cos(geo.beta1)-z^2);
    dfi1=atan(x/geo.Db);
    fprintf('\n x=%5.3f, y=%5.3f, z=%5.3f, dfi1=%5.3f, beta1=%5.3f',...
        x,y,z,dfi1*180/pi,geo.beta1*180/pi);
end

out=dfi1;

end

function dydx=streamlineode(t,z,C,geo)
tmp=vel(z(1)+1i*z(2),C,geo);
dydx=[tmp.u;tmp.v];
end

function [val,ter,dir]=streamlineevent(t,z,C,geo)

Ract=sqrt(z(1)^2+z(2)^2);
R_e=linspace(geo.Db,geo.D2,geo.N_r)/2;
val=Ract-R_e;
ter=zeros(size(val)); ter(end)=1;
dir=zeros(size(val));
end
