function data = optimize_ptx(U,E,alpha_set,eta,slotLength,R,noiseVar,capture,...
    SIC,AoI_thres,metric)

% function data = optimize_ptx(U,E,alpha_set,eta,slotLength,R,noise,capture,...
%    SIC,AoI_thres,metric)
%
% Optimize the transmission probabilities to minimize the average AoI, 
% minimize the AVP, or maximize the throughput of a slotted ALOHA system 
% with energy harvesting.
%
% The optimization parameter is ptx, a vector containing the probabilities
% that a device transmits a generated update if its battery level is [1 ... E]
%
% ﻿[1] Khac-Hoang Ngo, G. Durisi, A. Graell i Amat, A. Munari, and F. Lazaro, 
% “Age of information in slotted ALOHA with energy harvesting,” in IEEE 
% Global Communications Conference (Globecom), Kuala Lumpur, Malaysia, Dec. 
% 2023. [Online]. Available: https://research.chalmers.se/publication/537484/file/537484_Fulltext.pdf.
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
%   metric : select the metric, i.e., 'avgAoI', 'AVP', or 'throughput' to optimize
% 
% OUTPUTS: 
%   the structure 'data' that stores all system parameters and the
%   optimization results.

%% Debug mode?
debug = 0; 

if debug
    tic
    U = 30;   
    E = 3; 

    alpha_set = 1/U;      
    eta = 0.05;                
    
    slotLength = 100;     
    R = .8;                   
    noiseVar = db2pow(-20);
    
    capture = 0;
    SIC = 1;
    
    AoI_thres = 1000;

    metric = 'avgAoI'; %'AVP', 'throughput', or 'avgAoI'
end

%% Initialization
% Generate multiple initializations for the transmission probabilities.
% Here, we choose initial values of the form {0,1}^E, and also consider
% random initializations.
if E > 1
    ini_set = de2bi(0:2^E-1,E)';

    for ii = 1:100
        ini_set = [ini_set rand(E,1)];
    end
elseif E == 1
    warning('Let pi = 1. No need to optimize!')
    return
end
n_ini = size(ini_set,2);

%% Optimization
switch lower(metric)
    case 'avgaoi'
        metric_select = [1 0 0];
    case 'avp'
        metric_select = [0 1 0];
    case 'throughput'
        metric_select = [0 0 -1];
end

avgAoI  = zeros(length(alpha_set),1);
AVP     = zeros(length(alpha_set),1);
throughput = zeros(length(alpha_set),1);
ptx_opt = zeros(E,length(alpha_set));
w       = zeros(E+1,length(alpha_set));
pE_steady = zeros(E+1,length(alpha_set));

% optimize for each value of the update generation rate
for idxA = 1:length(alpha_set)
    alpha = alpha_set(idxA);

    fprintf('===========================\n');
    fprintf('U*Alpha = %1.1f  \n', alpha*U);
    fprintf('===========================\n');

    avgAoI_0    = zeros(n_ini,1);
    AVP_0       = zeros(n_ini,1);
    thru_0      = zeros(n_ini,1);
    ptx_opt_0   = zeros(E,n_ini);
    w_0         = zeros(E+1,n_ini);
    pE_steady_0 = zeros(E+1,n_ini);
    
    metric_eval = @(pii) approximation(U,E,alpha,eta,[0; pii],slotLength,R,...
                    noiseVar,capture,SIC,AoI_thres,metric_select);

    for idx = 1:n_ini
        % initialize the probabilities
        pi0 = ini_set(:,idx);
    
        % optimization for the selected metric
        [var_opt, obj_opt] = fminsearchcon(metric_eval,pi0,...
            zeros(size(pi0)),ones(size(pi0)));
    
        fprintf('----------\n');
        fprintf(metric)
        fprintf(':   %.4f\n', obj_opt)
        fprintf('optimal ptx: [0, ');
        fprintf('%1.4f, ', var_opt(1:end));
        fprintf('1]\n');
        
        % compute the other metrics
        [output,w_tmp,pE_steady_tmp] = approximation(U,E,alpha,eta,[0; var_opt],...
            slotLength,R,noiseVar,capture,SIC,AoI_thres,[1 1 1]);
        
        % save optimization results
        ptx_opt_0(:,idx) = var_opt;
        avgAoI_0(idx) = output(1);
        AVP_0(idx) = output(2);
        thru_0(idx) = output(3);
        w_0(:,idx) = w_tmp;
        pE_steady_0(:,idx) = pE_steady_tmp;
    end
    
    % choose the best optimization results over the initializations
    switch lower(metric)
        case 'avgaoi'
            [avgAoI_0,idxopt] = min(avgAoI_0);
            AVP_0 = AVP_0(idxopt);
            thru_0 = thru_0(idxopt);
            ptx_opt_0 = ptx_opt_0(:,idxopt);
            w_0 = w_0(:,idxopt);
            pE_steady_0 = pE_steady_0(:,idxopt);
        case 'avp'
            [AVP_0,idxopt] = min(AVP_0);
            avgAoI_0 = avgAoI_0(idxopt);
            thru_0 = thru_0(idxopt);
            ptx_opt_0 = ptx_opt_0(:,idxopt);
            w_0 = w_0(:,idxopt);
            pE_steady_0 = pE_steady_0(:,idxopt);
        case 'throughput'
            [thru_0,idxopt] = max(thru_0);
            AVP_0 = AVP_0(idxopt);
            avgAoI_0 = avgAoI_0(idxopt);
            ptx_opt_0 = ptx_opt_0(:,idxopt);
            w_0 = w_0(:,idxopt);
            pE_steady_0 = pE_steady_0(:,idxopt);
    end
    avgAoI(idxA) = avgAoI_0;
    AVP(idxA) = AVP_0;
    throughput(idxA) = thru_0;
    ptx_opt(:,idxA) = ptx_opt_0;
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
data.noiseVar = pow2db(noiseVar);
data.metric = metric;
data.avgAoI = avgAoI;
data.AVP = AVP;
data.throughput = throughput;
data.ptx_opt = ptx_opt;
data.w = w;
data.pE_steady = pE_steady;

if ~debug
    filename = ['optimize_ptx_'];
    if capture
        filename = [filename 'capture_'];
        if SIC
            filename = [filename 'SIC_'];
        end
    end
    filename = [filename metric '_U_' num2str(U) '_alphaU_' num2str(alpha_set*U) ...
        '_batCap_' num2str(E) '_eta_' num2str(eta) '_ageThres_' num2str(AoI_thres) ...
        '_noiseVar_' num2str(pow2db(noiseVar)) 'dB'];
    

    filename = [filename '.mat'];

    save(filename,'data','-v7.3');
else
    toc
    keyboard
end
end
