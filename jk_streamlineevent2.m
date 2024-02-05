function [val,ter,dir]=jk_streamlineevent2(t,z,C,S,geo)
val=z(1)^2+z(2)^2-(geo.D2/2)^2;
ter=1;
dir=0;
end