function [] = plotSamp(samp, param)
% Author: Rizwan Ahmad (ahmad.46@osu.edu)

p = param.p;
t = param.t;

% Plot the sampling pattern marginals
pp = find(samp==1);
[x,~]=ind2sub([p,t],pp);
x = diff(x(:));
x(x<0)=[];

figure; 
        subplot(211); plot((sum(samp,2))); xlabel('k_y'); ylabel('Number of samples');
        subplot(212); plot((sum(samp,1))); xlabel('frames'); ylabel('Number of samples');
%         subplot(313); plot(1:numel(x),x,1:numel(x), sort(x),'r'); xlabel('samples'); ylabel('k-space jumps');



% 2D Binary map
figure; imagesc(squeeze(samp(:,:))); colormap(gray); %axis('image');
title([param.typ, ', R_e = ' num2str(p*t/sum(samp(:)))]); ylabel('k_y'); xlabel('frame index');
% samp=logical(samp);


% 2D plot
% figure; 
% for i=1:param.t
%     k = (i-1)*param.tr + 1 : i*param.tr;
%     col = [rand,rand,rand];
%     col = col - 0.5*min(col);
%     col = (2/3 + 1/3*rand)*col/max(col);
%     plot(T(k),P(k),'o','markersize',param.sz,'color',[0 0 0],'markerfacecolor',col); axis('image'); hold on;
% end
% title([param.typ, ', Re = ' num2str(p*t/sum(samp(:)))]);
% axis([-floor(param.t/2), ceil(param.t/2)-1, -floor(param.p/2), floor(param.p/2)-1]);
