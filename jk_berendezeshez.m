rho=1000;
g=9.81;
H=1:0.1:40;
Qa=1:0.1:40;
Q=Qa/3600;
n_q = ones(length(Q),length(H));
n=2800;
kk=1;
for ii=1:length(Q)
    for jj=1:length(H)
        n_q(jj,ii)=n*Q(ii)^0.5/H(jj)^0.75;
        eta(jj,ii)=0.94-0.048*Q(ii)^(-0.32)-0.29*(log10(n_q(jj,ii)/44))^2;
        P_h(jj,ii)=Q(ii)*H(jj)*rho*g;
        P_o(jj,ii)=P_h(jj,ii)/eta(jj,ii);
        psi=(300/(270+n_q(jj,ii)))^(9/4);
        u_2=sqrt(2*g*H(jj)/psi);
        D_2(jj,ii)=u_2/(pi*n/60);
        c_2m=0.1011*sqrt(2*g*H(jj));
        b_2(jj,ii)=Q(ii)/(D_2(jj,ii)*pi*c_2m);
        if eta(jj,ii) > 0.5 && P_o(jj,ii) < 1400 && n_q(jj,ii) < 36;
            % figure(100)
            % plot(Q(ii)*3600,H(jj),'*');
            % hold on
        oD_2(kk)=D_2(jj,ii);
        ob_2(kk)=b_2(jj,ii); 
        oQ(kk)=Q(ii);
        oH(kk)=H(jj);
        nq(kk)=n_q(jj,ii);
        kk=kk+1;
        end
    end
end
% [oD_2max,idx1]=max(oD_2)
% oQ(idx1)*3600
% oH(idx1)
% nq(idx1)

% [ob_2max,idx2]=max(ob_2)
% oQ(idx2)*3600
% oH(idx2)
% oQmax=max(oQ)*3600
% min(oD_2)
% min(ob_2)
% ob_2(idx1)
% oD_2(idx2)

figure(1)
surf(Qa,H,n_q)
clim([35 36])

figure(2)
surf(Qa,H,eta)
colormap(flipud(parula))
clim([0.45 0.46])

figure(3)
surf(Qa,H,P_o)
clim([1380 1400])

figure(4)
surf(Qa,H,D_2)
clim([0.09 0.10])

figure(5)
surf(Qa,H,b_2)
clim([0.01 0.011])


%clim([0 20])


Q_c=5.7/3600

fi=0:60:360
A_A=fi*Q_c/360/(15.4769^-0.29)/((2*g*13.9)^0.5);
h_A=A_A/0.014
