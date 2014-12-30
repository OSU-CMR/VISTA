function [samp] = VISTA(param)
% Author: Rizwan Ahmad (ahmad.46@osu.edu)

typ   = param.typ;   % Type of sampling
p     = param.p;     % Number of phase encoding steps
t     = param.t;     % Number of frames
tf    = param.tf;    % Relative step-size in "time" dimension wrt to phase-encoding direction
W     = param.W;     % Scaling of time dimension
fs    = param.fs;    % To make time average fully sampled or not
uni   = param.uni;   % Iteration where sampling is reset to jittered uniform
R     = param.R;     % Net acceleration factor
nIter = param.nIter; % Number of iterations
alph  = param.alph;  % 0=uniform density, 0<alph<1 for variable density
sig   = param.sig;   % Std of the Gaussian envelope that defines variable density
s     = param.s;     % Exponent of force term
g     = param.g;     % At every g^th iteration, sampled are relocated to the nearest Cartesian grid
sz    = param.sz;    % Size of displayed samples
fl    = param.fl;    % At/beyond this iteration, start checking for fully-sampledness
ss    = param.ss;    % Gradient descent step-size
dsp   = param.dsp;   % Display the distribution at every dsp^th iteration
sd    = param.sd;    % Seed for random number generation

tr    = round(param.p/param.R); % Number of readout lines per frame (temporal resolution)
N     = tr*t; % Total number of samples after rounding off

%% Let's handle the special case of R = 1
if R == 1
    samp = ones(p, t);
    return;
end

%% Let's find UIS
if  strcmp(param.typ, 'UIS')
    %% Uniform sampling ===========================
    ptmp  = zeros(param.p,1); ptmp(round(1:param.R:param.p))=1;
    ttmp  = zeros(param.t,1); ttmp(round(1:param.R:param.t))=1;
    Top = toeplitz(ptmp,ttmp);
    ind = find(Top);
    ph  = rem(ind-1, param.p)-(floor(param.p/2));
    ti  = floor((ind-1)/param.p)-(floor(param.t)/2);
    [ph, ti] = dispdup(ph, ti, param);
    samp = zeros(p, t);
    ind  = round(p*(ti + floor(t/2)) + (ph + floor(p/2)+1));
    samp(ind) = 1;
    return
end

%% Use VRS as initialization for VISTA
p1=(-floor(p/2):ceil(p/2)-1)';
t1=[];
ind=0;
ti = zeros(tr*t,1);
ph = zeros(tr*t,1);
prob = (0.1 + alph/(1-alph+1e-10)*exp(-(p1).^2./(1*sig.^2)));

rng(sd);
tmpSd = round(1e6*rand(t,1)); % Seeds for random numbers
for i = -floor(t/2):ceil(t/2)-1
    a = find(t1==i);
    n_tmp = tr-numel(a);
    prob_tmp = prob;
    prob_tmp(a) = 0;
    p_tmp = randp(prob_tmp, tmpSd(i+floor(t/2)+1), n_tmp, 1) - floor(p/2)-1;
    ti(ind+1:ind+n_tmp) = i;
    ph(ind+1:ind+n_tmp) = p_tmp;
    ind = ind+n_tmp;
end
if strcmp(typ,'VRS')
    [ph, ti] = dispdup(ph, ti, param);
    samp = zeros(p, t);
    ind  = round(p*(ti + floor(t/2)) + (ph + floor(p/2)+1));
    samp(ind) = 1;
    return
end


% Displacement parameters
fprintf(['Computing VISTA, plese wait as it may take a while ...' '\n']);
tic,
stp = ones(1,nIter); % Gradient descent displacement
a = W*ones(1,nIter);  % Temporal axis scaling


f = 1+round(100*rand); % Figure index
dis_ext = zeros(N,1); % Extent of displacement
for i=1:nIter
    [ph,ti] = tile(ph(1:N), ti(1:N), param);
    for j=1:N
        
       
        % Distances -------------------------------------------------------
        dis= sqrt(abs(ph-ph(j)).^2 + abs(a(i)*(ti-ti(j))).^2);
        nanloc= dis==0;
        dis(nanloc)=inf;
        
        % Scaling ---------------------------------------------------------
        scl= 1 - alph*exp(-(ph).^2./(2*sig.^2)); 
        scl = scl + (1-scl(1));
        dscl = 1/sig^2 * alph.*ph(j).*exp(-(ph(j).^2)./(2*sig^2)); % Differentiation of scl wrt to "ph"

        % Force and resulting displacement --------------------------------
        fx = s*(     (ph(j)-ph) .* (scl(j)*scl./dis.^(s+2))) - dscl .* scl./dis.^s;
        fy = s*(a(i)^2*(ti(j)-ti) .* (scl(j)*scl./dis.^(s+2)))*tf;
        ph(j) = (ph(j)+ max(min(stp(i)*sum(fx), R/4),-R/4));
        ti(j) = (ti(j)+ max(min(stp(i)*sum(fy), R/4),-R/4));
                
        % Ensure that the samples stay in bounds --------------------------
        if ph(j) < -floor(p/2)-1/2 
            ph(j) = ph(j) + p;
        elseif ph(j)> (ceil(p/2)-1/2)
            ph(j) = ph(j) - p;
        end

        if ti(j) < -floor(t/2)-1/2 
            ti(j) = ti(j) + t;
        elseif ti(j)> (ceil(t/2)-1/2);
            ti(j) = ti(j) - t;
        end
        
        % Displacing samples to nearest Cartesian location
        if rem(i,g)==0 || i==nIter
            ph(j)=round(ph(j));
            ti(j)=round(ti(j));
        end
        
        % Measuing the displacement
        if i==2, dis_ext(j) = abs(stp(i)*sum(fx)); end
    end
    
    
    % Normalizing the step-size to a reasonable value
    if i==3
        stp = ss*(1+R/4)/median(dis_ext)*stp;
    end
    
    % At uni-th iteration, reset to jittered uniform sampling
    ti = ti(1:N);
    ph = ph(1:N);
    if i == uni
        tmp = zeros(tr, t);
        for k = 1:t
            tmp(:,k) = sort(ph((k-1)*tr+1:k*tr));
        end
        tmp = round(mean(tmp,2)); % Find average distances between adjacent phase-encoding samples
        ph = repmat(tmp, [t,1]); % Variable density sampling with "average" distances
        rng(sd);
        rndTmp = rand(t,1);
        for k=-floor(t/2):ceil(t/2)-1
            tmp = ti==k;
            ptmp = ph(tmp) + round(1/2*R^(1)*(rndTmp(k+floor(t/2)+1)-0.5)); % Add jitter
            ptmp(ptmp>ceil(p/2)-1) = ptmp(ptmp>ceil(p/2)-1)-p; % Ensure the samples don't run out of the k-t grid
            ptmp(ptmp<-floor(p/2)) = ptmp(ptmp<-floor(p/2))+p;
            ph(tmp) = ptmp;
        end
        % Temporarily stretch the time axis to avoid bad local minima
        a(i+1:end) = a(i+1:end) .* (1 + exp(-((i+1:nIter)-(i+1))/ceil(nIter/60))); 
        %figure; plot(a); 
    end
    
    % Displace the overlapping points so that one location has only one sample
    if rem(i,g)==0 || i==nIter;
        [ph, ti] = dispdup(ph(1:N), ti(1:N), param);
    end
    
    % Check/ensure time-average is fully sampled
    if (rem(i,g)==0 || i==nIter) && i>=fl 
        ph = ph(1:N);
        ti = ti(1:N);
        if fs==1; % Ensuring fully sampledness at average all
            [ph, ti] = fillK(ph, ti, ph, ti, param);
        elseif fs>1 % Ensuring fully sampledness for "fs" frames
            for m = 1:floor(t/fs)
                tmp = (m-1)*tr*fs + 1 : m*tr*fs;
                [ph, ti] = fillK(ph(tmp), ti(tmp), ph, ti, param);
            end
        end
    end    
    
    if (i==1 || rem(i,dsp)==0 || i==nIter) % When to diplay the distribution
        figure(f);
        plot(ti(1:N),ph(1:N),'s','markersize',sz,'color','red','markerfacecolor','red'); 
        axis('image');axis([-floor(t/2), ceil(t/2)-1, -floor(p/2), ceil(p/2)-1]);
        title(['Iter: ' num2str(i) ', Number of samples: ' num2str(N)]);
        hold off;
    end
end
[ph, ti] = dispdup(ph(1:N), ti(1:N), param);

% From indices to 2D binary mask
samp=zeros(p, t);
ind = round(p*(ti + floor(t/2)) + (ph + floor(p/2)+1));
samp(ind)=1;
fprintf('VISTA computed in %4.2f s\n', toc);



function [po,to] = tile(ph,ti,param)
% Replicate the sampling pattern in each direction. Probablity, this is 
% not an efficient way to impose periodic boundary condition because it 
% makes the problem size grow by 9 fold.

p = param.p; 
t = param.t; 
po=[ph; ph-p; ph; ph+p;     ph-p; ph+p;    ph-p; ph; ph+p];
to=[ti; ti-t; ti-t; ti-t;   ti; ti;        ti+t; ti+t; ti+t];
    
    


