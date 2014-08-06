function [ph, ti] = dispdup(ph, ti, param)
% If multiple samples occupy the same location, this routine displaces the 
% duplicate samples to the nearest vacant location so that there is no 
% more than one smaple per location on the k-t grid

% Author: Rizwan Ahmad (ahmad.46@osu.edu)

ph=ph+ceil((param.p+1)/2);
ti=ti+ceil((param.t+1)/2);
pt = (ti-1)*param.p + ph;

uniquept = unique(pt);
countOfpt = hist(pt,uniquept);
indexToRepeatedValue = (countOfpt~=1);
repeatedValues = uniquept(indexToRepeatedValue);
dupind=[];
for i=1:numel(repeatedValues)
    tmp = find(pt==repeatedValues(i));
    dupind = [dupind;tmp(2:end)]; % Indices of locations which have more than one sample
end

empind = (setdiff(1:param.p*param.t, pt))'; % Indices of vacant locations


for i=1:numel(dupind) % Let's go through all 'dupind' one by one
    newind = nearestvac(pt(dupind(i)),empind,param);
    pt(dupind(i))=newind;
    empind=setdiff(empind,newind);
end

[ph,ti] = ind2xy(pt,param.p);
ph=ph-ceil((param.p+1)/2);
ti=ti-ceil((param.t+1)/2);


% For a given 'dupind', this function finds the nearest vacant location 
% among 'empind'
function [n] = nearestvac(dupind, empind, param)
[x0, y0] = ind2xy(dupind, param.p);
[x ,y] = ind2xy(empind, param.p);
dis1 = (x-x0).^2;
dis2 = (y-y0).^2; dis2(dis2>eps)=inf;
dis = sqrt(dis1 + dis2);
[~,b] = min(dis);
n = empind(b);

% Index to (x,y)
function [x,y] = ind2xy(ind,X)
x = ind - floor((ind-1)/X)*X;
y = ceil(ind/X);

