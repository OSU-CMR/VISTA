function [Pacc, Tacc] = fillK(P, T, Pacc, Tacc, param)
% Ensures time-average of VISTA is fully sampled (except for the outer most
% region)
% Author: Rizwan Ahmad (ahmad.46@osu.edu)

p   = param.p;


% empty locations;
tmp = setdiff(-floor(p/2):ceil(p/2)-1,P);
[tmp2, ord] = sort(abs(tmp));
tmp2 = tmp2.*sign(tmp(ord)); % Sorted (from center-out) empty locations

while numel(tmp2)>0
    ind = tmp2(1); % the hole (PE location) to be filled
    
    can = P((sign(ind+eps)*P) > (sign(ind+eps)*ind)); % Find candidates which are on the boundary side of "ind"
    if isempty(can), break; end % If there is nothing on the boundary side
    Pcan = can(sign(ind+eps)*(can-ind) == min(sign(ind+eps)*(can-ind))); % P index of candidates that can be used to cover the void
    Tcan = T(P==Pcan(1));                                                % T index of candidates that can be used to cover the void

    U = inf;
    for i = 1:numel(Pcan)
        Ptmp = Pacc;
        Ttmp = Tacc;
        Ptmp(Pacc == Pcan(i) & Tacc == Tcan(i)) = ind; 
        Utmp = computeU(Ptmp,Ttmp,param); % Compute engergy U for the ith candidate
        if Utmp<=U
            slc = i;
            U = Utmp;
        end
    end

    Pacc(Pacc == Pcan(slc) & Tacc == Tcan(slc)) = ind; % Fill the hole with the appropriate candidate
    P(P == Pcan(slc) & T == Tcan(slc)) = ind; % Fill the hole with the approprate candidate
%     figure; plot(Tacc,Pacc,T,P+0.1,'--r');
    tmp = setdiff(-floor(p/2):ceil(p/2)-1,P);
    tmp = excludeOuter(tmp,p);
    [tmp2, ord] = sort(abs(tmp));
    tmp2 = tmp2.*sign(tmp(ord)); % Find new holes

% figure;
% plot(T,P,'s','markersize',sz,'color','red','markerfacecolor','red'); 
% axis('image');axis([-t/2, t/2-1, -p/2, p/2-1]);
    
end


function [tmp] = excludeOuter(tmp,p)
% Remove the "voids" that are contiguous to the boundary

tmp = sort(tmp,'descend');
cnd = abs(diff([ceil(p/2),tmp]));
if max(cnd)==1
    tmp = [];
else 
    tmp(1:find(cnd > 1,1)-1) = [];
end

tmp = sort(tmp,'ascend');
cnd = abs(diff([-floor(p/2)-1,tmp]));
if max(cnd)==1
    tmp = [];
else 
    tmp(1:find(cnd > 1,1)-1) = [];
end




