function data = optimize_ptx(U,E,alpha_set,eta,slotLength,R, ...
    noiseVar,capture,SIC,AoI_thres,metric)

% function data = optimize_ptx(U,E,alpha_set,eta,slotLength,R, ...
%    noise,capture,SIC,AoI_thres,metric)
%
% Optimize the transmission probabilities to minimize the average AoI, 
% minimize the AVP, or maximize the throughput of a slotted ALOHA system 
% with energy harvesting.
%
% The optimization parameter is ptx, a matrix whose (i,j)th entry is 
% the probability that a device transmits with j energy units if it has a 
% new update and its battery level is i \in [0,1,...,E].
%
% Reference:
%
% ﻿Khac-Hoang Ngo, G. Durisi, A. Munari, and F. Lazaro, and A. Graell i Amat,
% "Timely status updates in slotted ALOHA networks with energy harvesting," 
% submitted to IEEE Transactions on Communications, Apr. 2024.
%
% INPUTS: 
%   U     : number of devices
%   E     : battery capacity of a device
%   alpha_set : probability that a device has a new update in a slot
%   eta   : probability that a device harvests an energy unit in a slot
%   slotLength : number of channel uses in a slot
%   R     : transmission rate (bits/channel use)
%   noiseVar   : noise variance
%   capture    : indicator if decoding is with (1) or without (0) capture
%   SIC        : indicator if successive interference cancellation is
%                considered (1) or not (0)
%                Note: if capture = 0, SIC is not considered 
%   AoI_thres  : AoI threshold considered for AVP
%   metric : select the metric, i.e., 'avgAoI', 'AVP', or 'throughput', to optimize
% 
% OUTPUTS: 
%   the structure 'data' that stores all system parameters and the
%   optimization results.

addpath('./helpers')

%% Debug mode?
debug = 1; 
tic
if debug
    %% Parameters
    U = 1000;           % number of devices
    E = 8;              % battery capacity

    alpha_set = 1/U;    % update generation rate per device
    
    eta = 0.005;        % probability of harvesting an energy unit in a slot
    
    slotLength = 100;   % length of a slot
    R = .8;             % transmission rate
    
    AoI_thres = 10000;  % age threshold
    
    noiseVar = db2pow(-20);
    
    capture = 0;
    SIC = 1;

    metric = 'avgAoI';  % 'AVP', 'throughput', 'avgAoI'

end

warning off;

%% Initialization
% Generate multiple initializations for the transmission probabilities

minE_to_opt = 1;      % we assume that if the battery level is smaller than 
                      %   minE_to_opt, the device does not transmit. This
                      %   is to reduce the number of optimization variables

nvar = (E-minE_to_opt+1)*(E-minE_to_opt+2)/2;   % number of optimization variables

n_ini = 200;        % number of trials 
ini_set = [];       % initialization set
for idxini = 1:n_in
    ini_p = rand;
    for ii = 2:E-minE_to_opt+1
            ini_p = [ini_p; [randfixedsum(ii,1,rand,0,1)]];
    end
    ini_set = [ini_set ini_p];
end

condA = zeros(E-minE_to_opt+1,nvar);
tmp = 1;
for ii = 1:E-minE_to_opt+1
    condA(ii,tmp:tmp+ii-1) = 1;
    tmp = tmp + ii;
end
condb = ones(E-minE_to_opt+1,1);


%% Optimization
switch lower(metric)
    case 'avgaoi'
        metric_sel = [1 0 0];
    case 'avp'
        metric_sel = [0 1 0];
    case 'throughput'
        metric_sel = [0 0 -1];
end

avgAoI = zeros(length(alpha_set),1);
AVP = zeros(length(alpha_set),1);
thru = zeros(length(alpha_set),1);
ptx_opt = zeros(E+1,E+1,length(alpha_set));
w = zeros(E+1,length(alpha_set));
pE_steady = zeros(E+1,length(alpha_set));

% optimize for each value of the update generation rate
for idxA = 1:length(alpha_set)
    alpha = alpha_set(idxA);

    fprintf('===========================\n');
    fprintf('U*Alpha = %1.1f  \n', alpha*U);
    fprintf('===========================\n');

    avgAoI_0 = zeros(n_ini,1);
    AVP_0 = zeros(n_ini,1);
    thru_0 = zeros(n_ini,1);
    ptx_opt_0 = zeros(E+1,E+1,n_ini);
    w_0 = zeros(E+1,n_ini);
    pE_steady_0 = zeros(E+1,n_ini);
    
    metric_eval = @(pii) approximation(U,E,alpha,eta, ...
        optvar_converter(E,pii,minE_to_opt),[ones(E+1,1) zeros(E+1,E)], ...
        slotLength,R,noiseVar,capture,SIC,AoI_thres,metric_sel,1e5);
    
    valid_idx = zeros(n_ini,1);
    for idx = 1:n_ini   % use parfor to paralellize
        % initialize the probabilities
        pi0 = ini_set(:,idx);
    
        % optimization
        [var_opt, obj_opt] = fminsearchcon(metric_eval,pi0,...
            zeros(size(pi0)),ones(size(pi0)),condA,condb);
        
        if ~isempty(var_opt)
            ptx_opt_0(:,:,idx) = optvar_converter(E,var_opt,minE_to_opt);
        
            fprintf('----------\n');
            fprintf(metric)
            fprintf(':   %.4f\n', obj_opt)
            fprintf('optimal ptx: [%1.4f,',var_opt(1));
            fprintf('%1.4f, ', var_opt(2:end-1));
            fprintf('%1.4f]\n',var_opt(end));
            
            % compute the other metrics
            [output,w_tmp,pE_steady_tmp] = approximation(U,E,alpha,eta, ...
                ptx_opt_0(:,:,idx),[ones(E+1,1) zeros(E+1,E)],slotLength,R,...
                noiseVar,capture,SIC,AoI_thres,[1 1 1],1e5);
            avgAoI_0(idx) = output(1);
            AVP_0(idx) = output(2);
            thru_0(idx) = output(3);
            switch lower(metric)
                case 'avgaoi'
                    avgAoI_0(idx) = obj_opt;
                case 'avp'
                    AVP_0(idx) = obj_opt;
                case 'throughput'
                    thru_0(idx) = -obj_opt;
            end
            w_0(:,idx) = w_tmp;
            pE_steady_0(:,idx) = pE_steady_tmp;
            
            valid_idx(idx) = 1;
        end
    end

    % store only the valid results
    valid_idx = logical(valid_idx);
    avgAoI_0 = avgAoI_0(valid_idx);
    AVP_0 = AVP_0(valid_idx);
    thru_0 = thru_0(valid_idx);
    w_0 = w_0(:,valid_idx);
    pE_steady_0 = pE_steady_0(:,valid_idx);
    ptx_opt_0 = ptx_opt_0(:,:,valid_idx);

    % choose the best optimization results over the initializations
    switch lower(metric)
        case 'avgaoi'
            [avgAoI_0,idxopt] = min(avgAoI_0);
            AVP_0 = AVP_0(idxopt);
            thru_0 = thru_0(idxopt);
            ptx_opt_0 = ptx_opt_0(:,:,idxopt);
            w_0 = w_0(:,idxopt);
            pE_steady_0 = pE_steady_0(:,idxopt);
        case 'avp'
            [AVP_0,idxopt] = min(AVP_0);
            avgAoI_0 = avgAoI_0(idxopt);
            thru_0 = thru_0(idxopt);
            ptx_opt_0 = ptx_opt_0(:,:,idxopt);
            w_0 = w_0(:,idxopt);
            pE_steady_0 = pE_steady_0(:,idxopt);
        case 'throughput'
            [thru_0,idxopt] = max(thru_0);
            AVP_0 = AVP_0(idxopt);
            avgAoI_0 = avgAoI_0(idxopt);
            ptx_opt_0 = ptx_opt_0(:,:,idxopt);
            w_0 = w_0(:,idxopt);
            pE_steady_0 = pE_steady_0(:,idxopt);
    end

    avgAoI(idxA) = avgAoI_0;
    AVP(idxA) = AVP_0;
    thru(idxA) = thru_0;
    ptx_opt(:,:,idxA) = ptx_opt_0;
    w(:,idxA) = w_0;
    pE_steady(:,idxA) = pE_steady_0;
end

%% Save the results
data.simtime = toc;
data.U = U;
data.alpha = alpha_set;
data.E = E;
data.eta = eta;
data.slotLength = slotLength;
data.R = R;
data.AoI_thres = AoI_thres;
data.capture = capture;
data.SIC = SIC;
data.noiseVardB = pow2db(noiseVar);
data.metric = metric;
data.avgAoI = avgAoI;
data.AVP = AVP;
data.throughput = thru;
data.ptx_opt = ptx_opt;
data.w = w;
data.pE_steady = pE_steady;
data.nvar = nvar;
data.n_ini = n_ini;

if ~debug
    filename = ['optimize_ptx_'];
    if capture
        filename = [filename 'capture_'];
        if SIC
            filename = [filename 'SIC_'];
        end
    end
    filename = [filename metric '_U_' num2str(U) '_alphaU_' ...
        num2str(min(alpha_set)*U) 'to' num2str(max(alpha_set)*U)...
        '_E_' num2str(E) '_eta_' num2str(eta) '_ageThres_' num2str(AoI_thres) ...
        '_noiseVar_' num2str(pow2db(noiseVar)) 'dB'];
    

    filename = [filename '.mat'];

    save(filename,'data','-v7.3');
else
    toc
    % keyboard
end
end

%% Convert the optimization variable vector to matrix form
function var = optvar_converter(E,var_vec,minE_to_opt)
    var = zeros(E+1,E+1);
    ii = 1;
    for i = minE_to_opt+1:E+1
        for j = minE_to_opt+1:i
            var(i,j) = var_vec(ii);
            ii = ii + 1;
        end
    end
    var(:,1) = 1 - sum(var(:,2:end),2);
end