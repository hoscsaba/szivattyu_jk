function C=solve_for_C(geo,Cini,DO_PLOT)

options = optimoptions(@fsolve,'MaxFunctionEvaluations',10e3,'MaxIterations',15e3);

C=fsolve(@(x)fun_to_solve(x,geo,DO_PLOT),Cini);

%fini=fun_to_solve(Cini,geo,0);
%ff=fun_to_solve(C,geo,0);
%[Cini C]
%[norm(fini) norm(ff)]
%pause
end

function out=fun_to_solve(C,geo,DO_PLOT)

[H_act,~,veldata,geo]=jk_kompl_pot(C,geo,DO_PLOT);

for kk=1:length(geo.x_v)    
    z=geo.x_v(kk)+geo.y_v(kk)*1i;   
    tmp=vel(z,C,geo);
    out(kk,1)=dot([geo.n_x(kk) geo.n_y(kk)],[tmp.u;tmp.v]);
end
%out=out;

end

