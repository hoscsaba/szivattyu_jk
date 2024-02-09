function save_to_CFD(geo)
  fp=fopen("jk_A.txt","w");
  w=0.1;
  figure(1)
  for i=1:length(geo.x_g(:,1))
    r=sqrt(geo.x_g(i,1)^2+geo.y_g(i,1)^2);
    theta=atan(geo.y_g(i,1)/geo.x_g(i,1));
    if theta<0
      theta=theta+pi;
    end
    fprintf(fp,"\n %5.3e %5.3e %5.3e %5.3e",r*1000,theta,geo.b2*1000,w);
  plot(r,theta,'*'), hold on
  end
  fclose(fp);
end
