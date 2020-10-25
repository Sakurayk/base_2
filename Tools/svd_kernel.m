%Author: Masaki Onuki (masaki.o@msp-lab.org)
function h = svd_kernel(x,d)
n = length(x);
h = zeros(n,1);


for i=1:n
    
    h(i)=max(sqrt(x(i))-d,0)/sqrt(x(i));
    
end


h = h';
end
