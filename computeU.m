function U = computeU(P,T,param)
% Compute potential energy (U) of the distribution
% Author: Rizwan Ahmad (ahmad.46@osu.edu)

N    = numel(P);
alph = param.alph;
sig  = param.sig;
a    = param.W;
s    = param.s;

U = 0;
for i = 1:N
%     for j = [1:i-1, i+1:N]
%         U = U + 1/2*((1 - alph*exp(-(P(i)).^2./(2*sig.^2))) * (1 - alph*exp(-(P(j)).^2./(2*sig.^2)))) ./...
%             ((P(i)-P(j))^2 + (a*(T(i)-T(j)))^2)^(e/2);
%     end
    
    k = [1:i-1, i+1:N];
        U = U + sum(1/2*((1 - alph*exp(-(P(i)).^2./(2*sig.^2))) .* (1 - alph*exp(-(P(k)).^2./(2*sig.^2)))) ./...
            ((P(i)-P(k)).^2 + (a*(T(i)-T(k))).^2).^(s/2));
    
end
