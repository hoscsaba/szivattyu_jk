function geo=jk_postprocess(geo)
    figure(123),clf, hold on
    subplot(2,3,[1,2,4,5])

    plot(geo.x_g,geo.y_g,'k','LineWidth',1)
    hold on
    RR=1.01;
    mul=1/200;
    r_b=geo.Db/2;
    r_k=geo.D2/2;
    phi_vec=linspace(0,2*pi,7);

    fi_ini=pi*linspace(-30,10,20)/180;
    ode_options=odeset('Events',@(t,z) jk_streamlineevent(t,z,geo.C,geo.S,geo));
    for ii=1:length(fi_ini)
        xini=0.8*r_b*cos(fi_ini(ii));
        yini=0.8*r_b*sin(fi_ini(ii));
        [ts,xsys]=ode45(@(t,z) jk_streamlineode(t,z,geo.C,geo.S,geo),...
            [0 0.1],[xini;yini],ode_options); %DE megoldó

        plot(xsys(:,1),xsys(:,2),'r');
    end

    hold off

    title(['N_r=',num2str(geo.N_r),', Q=',num2str(round(geo.QQ*3600)),...
        ' m3/h, H=',num2str(round(geo.HH*100)/100),'m'])
    xylim=r_k*2;
    xlim([-xylim,xylim]), ylim([-xylim,xylim]);
    axis("equal");

    drawnow
    
end