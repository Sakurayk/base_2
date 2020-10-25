%Author: Masaki Onuki (masaki.o@msp-lab.org)
function r=chebyshev_oprator(I,B,c,Range)
M = numel(c);
assert(all(M>=2));

% a1=Range(2)/2;

B_hat = (2/Range(2))*B - I;
CP_old=I; %Initialized Chebyshev polynomial (k=0 in)
CP_cur=B_hat; % j=1;

r = .5*c(1)*CP_old + c(2)*CP_cur;
for k=2:M
    CP_new = 2*B_hat*CP_cur-CP_old;
    if 1+k<=M
        r=r+c(k+1)*CP_new;
    end
    CP_old=CP_cur;
    CP_cur=CP_new;
    clear CP_new
end