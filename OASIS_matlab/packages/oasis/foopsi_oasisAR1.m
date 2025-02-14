function [c, s, b, g, active_set] = foopsi_oasisAR1(y, g, lam, smin, optimize_b,...
    optimize_g, decimate, maxIter, tau_range, gmax)
%% Infer the most likely discretized spike train underlying an AR(1) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g c_{t-1} >= 0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
%impulse response.
%   lam:  scalar, sparsity penalty parameter
%   optimize_b: bool, optimize baseline if True
%   optimize_g: integer, number of large, isolated events to consider for
%       optimizing g
%   decimate: int, decimation factor for estimating hyper-parameters faster
%       on decimated data
%   maxIter:  int, maximum number of iterations
%   active_set: npool x 4 matrix, warm stared active sets
%   gmax:  scalar, maximum value of g
%   tau_range: [tau_min, tau_max]
%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)
%   b: scalar, fluorescence baseline
%   g: scalar, parameter of the AR(1) process
%   active_set: npool x 4 matrix, active sets

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging


%% input arguments
y = reshape(y, [], 1);
T = length(y);

if ~exist('g', 'var') || isempty(g)
    g = estimate_time_constant(y, 1);
end
if ~exist('lam', 'var') || isempty(lam);   lam = 0; end
if ~exist('smin', 'var') || isempty(smin)
    smin = 0;
elseif smin<0
    smin = abs(smin) * GetSn(y);     % use a threshold that is proportional to the noise level
end
if ~exist('optimize_b', 'var') || isempty(optimize_b)
    optimize_b = false;
end
if ~exist('optimize_g', 'var') || isempty(optimize_g)
    optimize_g = 0;
end
if ~exist('decimate', 'var') || isempty(decimate)
    decimate = 1;
else
    decimate = max(1, round(decimate));
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 10;
end
if ~exist('tau_range', 'var') || isempty(tau_range)
   g_range = [0, 1];
else
   g_range = exp(-1./tau_range); 
   g = min(max(g, g_range(1)), g_range(2)); 
end

% change parameters due to downsampling
if decimate>1
    decimate = 1;  %#ok<NASGU>
    disp('to be done');
    %     fluo = y;
    %     y = resample(y, 1, decimate);
    %     g = g^decimate;
    %     thresh = thresh / decimate / decimate;
    %     T = length(y);
end

%% optimize parameters
if ~optimize_b   %% don't optimize the baseline b
    %% initialization
    b = 0;
    [solution, spks, active_set] = oasisAR1(y, g, lam, smin);
    
    %% iteratively update parameters g
    if  optimize_g     % update g
        [solution, active_set, g, spks] = update_g(y, active_set,lam, smin, g_range);
    end
else
    %% initialization
    b = quantile(y, 0.15);
    [solution, spks, active_set] = oasisAR1(y-b, g, lam, smin);
    
    %% iteratively update parameters g and b
    for m=1:maxIter
        b = mean(y-solution);
        if  optimize_g     % update g
            if size(active_set, 1) == 0
                break;
            end
            g0 = g;
            if g>gmax  % spike counts are too small. stop
                g = estimate_time_constant(y, 1);
                [solution, spks, active_set] = oasisAR1(y-b, g, lam, smin);
                break;
            end
            [solution, active_set, g, spks] = update_g(y-b, active_set,lam, smin, g_range);
            if abs(g-g0)/g0 < 1e-3  % g is converged
                optimize_g = false;
            end
        else
            break;
        end
    end
end
c = solution;
s = spks;
end
%update the AR coefficient: g
function [c, active_set, g, s] = update_g(y, active_set, lam, smin, g_range)
%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   active_set: npools*4 matrix, previous active sets.
%   lam:  scalar, curret value of sparsity penalty parameter lambda.

%% outputs
%   c: T*1 vector
%   s: T*1 vector, spike train
%   active_set: npool x 4 matrix, active sets
%   g: scalar

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

%% initialization
if ~exist('g_range', 'var') || isempty(g_range)
   g_range = [0, 1];  
end
len_active_set = size(active_set, 1);  %number of active sets
y = reshape(y,[],1);    % fluorescence data
maxl = max(active_set(:, 4));   % maximum ISI
c = zeros(size(y));     % the optimal denoised trace

%% find the optimal g and get the warm started active_set
g = fminbnd(@rss_g, g_range(1), g_range(2));
yp = y - lam*(1-g);
for m=1:len_active_set
    tmp_h = exp(log(g)*(0:maxl)');   % response kernel
    tmp_hh = cumsum(tmp_h.*tmp_h);        % hh(k) = h(1:k)'*h(1:k)
    li = active_set(m, 4);
    ti = active_set(m, 3);
    idx = ti:(ti+li-1);
    active_set(m,1) = (yp(idx))'*tmp_h(1:li);
    active_set(m,2) = tmp_hh(li);
end
[c,s,active_set] = oasisAR1(y, g, lam, smin, active_set);

%% nested functions
    function rss = rss_g(g)
        h = exp(log(g)*(0:maxl)');   % response kernel
        hh = cumsum(h.*h);        % hh(k) = h(1:k)'*h(1:k)
        yp = y - lam*(1-g);     % include the penalty term
        for ii=1:len_active_set
            li = active_set(ii, 4);
            ti = active_set(ii, 3);
            idx = ti:(ti+li-1);
            tmp_v = max(yp(idx)' * h(1:li) / hh(li), 0);
            c(idx) = tmp_v*h(1:li);
        end
        res = y-c;
        rss = res'*res;     % residual sum of squares
    end
end