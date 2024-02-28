function [QQ,HH,veldata,geo]=jk_main_get_QH(C,S,geo,DO_PLOT)

[QQ,HH,veldata.phi_ker_vec,veldata.c_k_u_vec,veldata.c_k_m_vec]=Get_QH(C,S,geo);

veldata.Gamma_lapat=Get_Circulation(geo,C,S);

end

function out=Get_Circulation(geo,C,S)

dx=1e-2;

ds=geo.t_arclength(end)/(geo.N_r-1);
%fprintf('\n');
for ii=1:geo.N_r-1
    z=geo.x_v(ii)+1i*geo.y_v(ii);
    vf=jk_vel(z + (geo.n_x(ii)+1i*geo.n_y(ii))*dx,C,S,geo);
    va=jk_vel(z - (geo.n_x(ii)+1i*geo.n_y(ii))*dx,C,S,geo);
    %fprintf('\n ii=%g, ds=%5.3f, vf=%5.3e, va=%5.3e, dv=%5.3e',...
    %   ii,ds,sqrt(vf.u^2+vf.v^2),sqrt(va.u^2+va.v^2),...
    %    sqrt(vf.u^2+vf.v^2)-sqrt(va.u^2+va.v^2));
    gg(ii)=ds*(sqrt(vf.u^2+vf.v^2)-sqrt(va.u^2+va.v^2));    
end
out=ds*sum(gg);

end


function [Q,H,phi_ker_vec,c_k_u_vec,c_k_m_vec]=Get_QH(C,S,geo)

RR=1.01;
r_k=geo.D2/2;
phi_ker_vec=linspace(0,2*pi,50*geo.N_lapat);
for i_phi=1:length(phi_ker_vec)
    phi_ker=phi_ker_vec(i_phi);
    n_k=[cos(phi_ker);sin(phi_ker)];
    t_k=[sin(phi_ker);-cos(phi_ker)];
    x_r_k=RR*r_k*cos(phi_ker);
    y_r_k=RR*r_k*sin(phi_ker);
    velo=jk_vel(x_r_k+y_r_k*1i,C,S,geo);
    w_k=[velo.u;velo.v];
    u_k=RR*r_k*geo.omega*t_k;
    c_k=w_k+u_k;
    c_k_u_vec(i_phi)=dot(c_k,t_k);
    c_k_m_vec(i_phi)=dot(c_k,n_k);
end

Q=mean(c_k_m_vec)*2*r_k*pi*geo.b2;
H=mean(c_k_u_vec)*norm(u_k)/9.81;
end
