%Author: Masaki Onuki (masaki.o@msp-lab.org)
function c=chebyshev_coefficient(h, Approx_order,Range)

a1=(Range(2)-Range(1))/2;
a2=(Range(2)+Range(1))/2;

N = Approx_order+1;

for j=1:Approx_order+1
    c(j)=sum (h(a1* cos( (pi*((1:N)-0.5))/N) + a2).*cos(pi*(j-1)*((1:N)-.5)/N) )*2/N;
end
