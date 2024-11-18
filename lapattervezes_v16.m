function lapattervezes_v16
clear all, close all

%% Main inputs
Q_target=2.83333/1000; H_target=16;
n=2345; g=9.81;
N_lapat=8; N_r=40; % min 4!!!
fname_prefix='jk_1';

%% Calculating the main parameters of the pump
nq=n*Q_target^0.5/H_target^0.75; 
psi=(300/(270+nq))^(9/4);
%u2=sqrt(2*g*H_target/psi);
%u2=sqrt(2*g*H_target/psi)*1.011;
%D2=u2/(pi*(n/60));
D2=0.138;
%u2=D2*pi*n/60;
Db=D2*0.4;
omega=2*pi*n/60;
u1=Db*pi*n/60; u2=D2*pi*n/60;
geo.u2=u2;
epszilon=0.0188*nq^(2/3); c1=epszilon*sqrt(2*g*H_target);
beta1=atan(c1/u1);
k2m=0.06+0.00195*nq;
c2u=H_target*g/u2; c2m=k2m*sqrt(2*g*H_target);
geo.beta2=atan(c2m/(u2-c2u));
%b2=1.05*Q_target/(D2*pi*c2m)
b2=0.0043;
%b1=Q_target/(Db*pi*c1);



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

%% Run a series of computations
%Percentages of H/Q
iterv_A=0;           %A_C initialization
iterv_i=1;              %running index
while iterv_A<6
A_C=iterv_A/(N_r);
A_S=0.0033e-4;
[ff,geo]=obj(A_C,A_S,geo,0,iterv_A);

fname=[fname_prefix,'_Q_',num2str(round(Q_target*3600)),'m3ph_H_',...
     num2str(round(geo.H_target)),'m.mat']
 save(fname,"geo");
 iterv_diff(iterv_i)=(geo.HH-H_target)^2; %difference square
 iterv_Av(iterv_i)=iterv_A;
 if (iterv_i>1)
     iterv_diff(iterv_i);
     iterv_diff(iterv_i-1);
     if (iterv_diff(iterv_i)>iterv_diff(iterv_i-1))
         %A_C=iterv_A-0.1/(N_r);
         break
     end
 end
% geo.N_r-1
iterv_A=iterv_A+0.01;
iterv_i=iterv_i+1;
end
A_C=iterv_Av(end-1)/(N_r);
iterv_Av(end-1);
[ff,geo]=obj(A_C,A_S,geo,1,iterv_A);
%%
%save_to_CFX(geo)
jk_postprocess(geo);

%%

% figure(153)
% plot(iterv_Av,iterv_diff,'LineWidth',2)
% hold on
% plot(iterv_Av(end-1),iterv_diff(end-1),'*','MarkerSize',12,'LineWidth',1.5, ...
%     'Color', [0.8500 0.3250 0.0980])
% for iterii=2:15:62
% plot(iterv_Av(iterii),iterv_diff(iterii),'*','MarkerSize',12,'LineWidth',1.5, ...
%     'Color', [0.4660 0.6740 0.1880])
% end
% xlabel('A_g') 
% ylabel('(H-H_{terv})^2 (m^2)') 


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
    if xi(i)<0.5
        C(i)=sin(acos(1-2*xi(i)));
    else
        C(i)=sin(acos(2*xi(i)-1));
    end
    if i==1
    Int_C=Int_C+C(i)*(xi(i));
    else
    Int_C=Int_C+C(i)*(xi(i)-xi(i-1));
    end
end
% figure(109)
% subplot(2,3,6)
% plot (xi,C,'-o','LineWidth',2)
% ylabel('\Gamma_0') 
% xlabel('\xi') 
%Gamma_r=Int_C*geo.t_arclength(end)*A_C/geo.b2
%Gamma_r/geo.Gamma_lapat_elm
%A=9.81*geo.H_target*2*pi*geo.b2/geo.N_lapat/geo.omega/geo.t_arclength(end)/Int_C
C=A_C*C;
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
function [out,geo]=obj(A_C,A_S,geo,DO_PLOT,iterv_A)
xi=geo.loc_c/geo.t_arclength(end);

delta_d_phi=1e5;
iter=1;
ITER_MAX=50;
iter_k=1;
while (delta_d_phi>0.01) && (iter<ITER_MAX)
    C=get_C(A_C,xi,geo);
    S=get_S(A_S,xi,geo);
    d_phi_old=geo.d_phi;
    alpha=2*pi/geo.N_lapat/2;
    tmax=0.05;
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

                figure(105) 
                plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
                axis('equal')

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

        [QQ,HH,veldata,geo]=jk_main_get_QH(C,S,geo,DO_PLOT);
        geo.QQ=QQ;
        geo.HH=HH;
        geo.veldata=veldata;

        err1=(geo.H_target-HH)^2;
        err2=0*std(geo.veldata.c_k_u_vec);
    end
    fprintf('\n\t iter #%d, norm(d_phi change) = %5.3e',iter,delta_d_phi);

    iter=iter+1;

    figure(105)

    plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
    title(['Q=',num2str(round(QQ*60000)),...
        'l/min, H=',num2str(round(10*HH)/10),'m, iter #',num2str(iter)])
    geo.d_phi';
    axis('equal')

    % if (iterv_A==0 && rem(iter_k,2)==1)
    %     figure(107)
    %     subplot(2,3,(iter_k+1)/2)
    %     plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
    %     title(iter_k+". iteráció")
    %     axis('equal')
    % else if (iterv_A==0 && iter_k==10)
    %     figure(107)
    %     subplot(2,3,(iter_k+2)/2)
    %     plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
    %     title(iter_k+". iteráció")
    %     axis('equal')
    % end
    % end

    % if (round(rem(iterv_A*100,15))==2)
    %     iterv_A*100
    %     (iterv_A*100-2)/15+1
    %     figure(108)
    %     subplot(2,3,round((iterv_A*100-2)/15+1))
    %     plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
    %     title(iterv_A*100+". lépés, A_g="+iterv_A)
    %     axis('equal')
    % end
    % if (round(iterv_A*100)==66)
    %     figure(108)
    %     subplot(2,3,6)
    %     plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
    %     title(iterv_A*100+". lépés, A_g="+iterv_A)
    %     axis('equal')
    % end



iter_k=iter_k+1;
%pause
if iter==ITER_MAX
    warning("Blade iteration not converging!");
    pause
end
end

out=0*err1+5*err2;
geo.C=get_C(A_C,xi,geo);
geo.S=get_S(A_S,xi,geo);
if (DO_PLOT==1)
        figure(100)
        subplot(2,3,[1,2,4,5])
        plot(geo.x_g, geo.y_g,'k',xsys(:,1),xsys(:,2),'r')
        subplot(2,3,3)
        xx=linspace(0,1);
        plot(xx,get_C(A_C,xx,geo))
        subplot(2,3,6)
        plot(xx,get_S(A_S,xx,geo))
end

        geo_matrix=[geo.x_g(:,1),geo.y_g(:,1)];
        writematrix(geo_matrix,'geo_matrix.xlsx','Sheet',1)


%% Send the geometry to postprocessing
%%ITT
%geo=jk_postprocess(geo);
%save_to_CFX(geo)
%geo.t_arclength(end)
end




function [val,ter,dir]=streamlineevent(t,z,C,S,geo)

Ract=sqrt(z(1)^2+z(2)^2);
R_e=linspace(geo.Db,geo.D2,geo.N_r)/2;
val=Ract-R_e;
ter=zeros(size(val)); ter(end)=1;
dir=zeros(size(val));
end
