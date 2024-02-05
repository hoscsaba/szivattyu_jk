xt=6;
x=0:0.1:xt;
xv=2;
xvv=xt-xv;
y=zeros(size(x));
for ii=1:length(x);
y(ii)=xnegyzet(x(ii),xv,xvv,xt);
end
plot(x,y);

fun = @(x) exp(-x.^2).*log(x).^2;
q = integral(fun,0,Inf)

function y=xnegyzet(x,xv,xvv,xt);
if x<xv
    y=sqrt(x)*(xv-x)^2;
elseif x>xvv
    y=-sqrt(-x+xt)*(xv-xt+x)^2;
else
    y=0;
end
end

