function geo=jk_build_geo2(geo)

N_r=length(geo.d_phi)+1;

phi(1)=0;
for ii=1:N_r-1
    phi(ii+1)=phi(ii)+geo.d_phi(ii);
end

r_b=geo.Db/2;
r_k=geo.D2/2;

% x_g,y_g -> geometria
for jj=1:geo.N_lapat
    for ii=1:N_r
        geo.x_g(ii,jj)=(r_b+(r_k-r_b)*(ii-1)/(N_r-1))*cos(phi(ii)+(jj-1)*2*pi/geo.N_lapat);
        geo.y_g(ii,jj)=(r_b+(r_k-r_b)*(ii-1)/(N_r-1))*sin(phi(ii)+(jj-1)*2*pi/geo.N_lapat);
    end
end

geo.t_arclength(1)=0;
for ii=1:N_r-1
    geo.t_arclength(ii+1)=geo.t_arclength(ii)+...
        norm([geo.x_g(ii+1,1)-geo.x_g(ii,1) geo.y_g(ii+1,1)-geo.y_g(ii,1)]);
end

% x_c,y_c -> cirkulaciok helye
for jj=1:geo.N_lapat
    for ii=1:geo.N_r-1
        geo.x_c(ii,jj)=(geo.x_g(ii,jj)+geo.x_g(ii+1,jj))/2;
        geo.y_c(ii,jj)=(geo.y_g(ii,jj)+geo.y_g(ii+1,jj))/2;
        if ii==1
            geo.loc_c(1)=norm([geo.x_c(1,1)-geo.x_g(1,1), geo.y_c(1,1)-geo.y_g(1,1)]);
        else
            geo.loc_c(ii)=geo.loc_c(ii-1)+...
                norm([geo.x_c(ii,1)-geo.x_c(ii-1,1) geo.y_c(ii,1)-geo.y_c(ii-1,1)]);
        end
    end
end

% x_v,y_v -> sebessegek kiertekelesenek a helyei
% tmp=linspace(0,geo.t_arclength(end),geo.N_r-1);

geo.x_v=geo.x_g(:,1);%(2:end,1);
geo.y_v=geo.y_g(:,1);%(2:end,1);

for ii=1:length(geo.x_v)
    if ii==1
        x1_t=[geo.x_v(1,1) geo.y_v(1,1)];
        x2_t=[geo.x_v(2,1) geo.y_v(2,1)];
    elseif ii==N_r
        x1_t=[geo.x_v(end-1,1) geo.y_v(end-1,1)];
        x2_t=[geo.x_v(end,1) geo.y_v(end,1)];
    else
        x1_t=[geo.x_v(ii-1) geo.y_v(ii-1)];
        x2_t=[geo.x_v(ii+1) geo.y_v(ii+1)];
    end
    t_x(ii)=x2_t(1)-x1_t(1);
    t_y(ii)=x2_t(2)-x1_t(2);
    t_norm=norm([t_x(ii) t_y(ii)]);
    geo.t_x(ii)=t_x(ii)/t_norm;
    geo.t_y(ii)=t_y(ii)/t_norm;
    geo.n_x(ii)=t_y(ii);
    geo.n_y(ii)=-t_x(ii);
end
geo.x_v=geo.x_v(2:end);
geo.y_v=geo.y_v(2:end);
geo.t_x=geo.t_x(2:end);
geo.t_y=geo.t_y(2:end);
geo.n_x=geo.n_x(2:end);
geo.n_y=geo.n_y(2:end);

end
