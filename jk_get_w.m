function [s,r,w,p,p_cp]=jk_get_w(xy,geo)
s(1)=0; %s means actual arclength of the blade
for i=1:length(xy(:,1))
    z=xy(i,1)+1i*xy(i,2);
    if i>1
        dz=[xy(i,1)-xy(i-1,1), xy(i,2)-xy(i-1,2)];
        s(i)=s(i-1)+norm(dz);
    end
    r(i)=norm(z);
    tmp=jk_vel(z,geo.C,geo.S,geo);
    wx=tmp.u; wy=tmp.v;
    w_v=[wx;wy];
    w(i)=norm(w_v);
    fi=angle(z); n=[-sin(fi); cos(fi)];
    u_v=norm(z)*geo.omega*n;
    u(i)=norm(u_v);
    c(i)=norm(w_v+u_v);

    if i==1
        w1=w(i); u1=u(i);
    end
    p(i)=1000 * ( (w1^2-w(i)^2)/2-(u1^2-u(i)^2)/2);

    tx=tmp.u;
    ty=tmp.v;
    n_v=[-ty;tx]/norm([tx ty]);

    % F_cp=m*r*omega^2, m=h*b2*ds*rho, p_cp=F_cp/(b2*ds)=rho*h*r*omega^2
    t_z_x=real(z);
    t_z_y=imag(z);
    v_cp=[t_z_x;t_z_y]/norm(z);
    p_cp_v=7860*0.05*r(i)*(geo.omega)^2*v_cp;

    p_cp(i)=dot(p_cp_v,n_v);
    w

    %fprintf('\n D/D1=%5.3f, c=%5.3f, u=%5.3f, p=%5.3f vom',r(i)/(geo.Db/2),norm(c),norm(u_v),p(i)/9.81/geo.rho);
end