function [output, w, pE_steady] = approximation(U,E,alpha,eta,ptx,slotLength,R,...
                    noiseVar,capture,SIC,AoI_thres,metric_select)
% function output = approximation(U,E,alpha,eta,ptx,slotLength,R,...
%                    noiseVar,capture,SIC,AoI_thres,metric_select,nMC)
%
% Evaluate approximations of the average age of information (AoI), age ...
% violation probability (AVP), and throughput of a slotted ALOHA system 
% with energy harvesting.
%
% ﻿[1] Khac-Hoang Ngo, G. Durisi, A. Graell i Amat, A. Munari, and F. Lazaro, 
% “Age of information in slotted ALOHA with energy harvesting,” in IEEE 
% Global Communications Conference (Globecom), Kuala Lumpur, Malaysia, Dec. 
% 2023. [Online]. Available: https://research.chalmers.se/publication/537484/file/537484_Fulltext.pdf.
%
% INPUTS: 
%   U     : number of devices
%   E     : battery capacity of a device
%   alpha : probability that a device has a new update in a slot
%   eta   : probability that a device harvests an energy unit in a slot
%   ptx   : a vector contains the probabilities that a device transmits a
%           generated update if its battery level is [0 1 ... E]
%   slotLength : number of channel uses in a slot
%   R     : transmission rate (bits/channel use)
%   noiseVar   : noise variance
%   capture    : indicator if decoding is with (1) or without (0) capture
%   SIC        : indicator if successive interference cancellation is
%                considered (1) or not (0)
%                Note: if capture = 0, SIC is not considered 
%   AoI_thres  : AoI threshold considered for AVP
%   metric_select : select the metric to return, set to [1 1 1] to return
%               all 3 metrics
%   nMC     : number of Monte-Carlo realizations for the power profiles of the devices 
% 
% OUTPUTS: output = [avgAoI; AVP; throughput], where
%   avgAoI : average AoI
%   AVP    : age violation probability
%   throughput : average number of packet delivered per slot
%   NOTE: some of the three metrics are not returned if the corresponding 
%         entry of metric_select is set to 0 

%% Debug mode?
debug = 0;

if debug
    tic
    U = 1000;   
    E = 8;  

    alpha = 2/U;     
    eta = 0.005;                 
    
    ptx = zeros(E+1,1);         
    ptx(2:E+1) = 1;    
    
    slotLength = 100;      
    R = .8;               
    noiseVar = db2pow(-20);

    capture = 0;
    SIC = 1;

    AoI_thres = 10000;  

    metric_select = [1 1 1];
end

output = [];

% probability that a device accesses the channel (i.e., generates a new 
% update and transmits it) 
rho = alpha*ptx;

%% Battery evolution
% --- Transition probabilities between battery levels, pE_transition
pE_transition = zeros(E+1,E+1);
for i = 0:E
for j = 0:E
    if i==j 
        if i  < E, pE_transition(i+1,j+1) = (1-eta)*(1-rho(i+1)); end
        if i == E, pE_transition(i+1,j+1) = (1-rho(i+1)); end
    elseif j == i + 1
        pE_transition(i+1,j+1) = eta*(1-rho(i+1)); 
    elseif j== 0
        pE_transition(i+1,j+1) = rho(i+1); 
    end
end
end

% --- Steady-state PMF of the battery level, pE_steady
Qhat = pE_transition';
Qhat(1,:) = 1;
for b = 2:E+1
    Qhat(b,b) = Qhat(b,b) - 1;
end

Qhat_inv = pinv(Qhat);
pE_steady = max(Qhat_inv(:,1),0); 
pE_steady = pE_steady/sum(pE_steady); % this is to avoid nummerical issue

%% Successful delivery probability
% --- Function to compute the error probability in AWGN channel
if capture
    dispersion = @(P,Pint,Pint2) log2(exp(1))^2.*(P.^2.*(1+2.*Pint + Pint.^2 ...
                - Pint2) + 2*P.*(Pint+1).^3)./(2*(Pint + 1).^2.*(P+Pint+1).^2);
    err = @(P,Pint,Pint2) qfunc(sqrt(slotLength./dispersion(P,Pint,Pint2))...
            .*(log2(1+P./(Pint+1))/2 - R));
    % precompute the error prob. in case of no interference (to speed up
    % the computation)
    err_precomp = err((1:E)/slotLength/noiseVar,0,0);
else
    dispersion = @(SNR) log2(exp(1))^2.*SNR.*(SNR+2)./(1+SNR).^2/2;
    err = @(SNR) qfunc(sqrt(slotLength./dispersion(SNR)).*(log2(1+SNR)/2 - R));
end

% --- Compute the successful delivery probability corresponding to each 
% battery level, which involves an expectation over a multinomial distribution

% initialize the successfull delivery probability 
w = zeros(E+1,1);

% probability that a device transmits with a given power
rhotmp = diag(rho); rhotmp(:,1) = 1 - sum(rhotmp(:,2:end),2);
pE_transmit = pE_steady'*rhotmp;

% number of Monte-Carlo simulation for the power profile of all devices
nMC = 1e5;

% now, compute the successfull delivery probability, w
if capture   
    if SIC % decoding with capture and SIC
        % generate realizations of the power profile of the active devices 
        Ue_trans = mnrnd(U,pE_transmit,nMC);
        Ue_trans(:,1) = [];

        % compute the success probability for each realization
        w_SIC_tmp = zeros(E,1);
        for idxMC = 1:nMC
            if sum(Ue_trans(idxMC,:))
                % transmit powers of each active devices
                powers = flip(repelem(1:E,Ue_trans(idxMC,:)));

                % compute the success probability of each active device
                w_tmp = compute_w_capture(powers,slotLength,noiseVar,1,err,err_precomp);
                for ii = 1:E
                    w_SIC_tmp(ii) = w_SIC_tmp(ii) + sum(w_tmp(powers == ii));
                end
            end
        end

        % average over the active devices
        w(2:end) = w_SIC_tmp./(sum(Ue_trans,1)');
        w(isnan(w)) = 0;
    else % decoding with capture and without SIC
        % transmit power profile of the other devices
        Ue_trans = mnrnd(U-1,pE_transmit,nMC);
        Ue_trans(:,1) = [];

        % successful delivery probabilities
        P = [1:E]/slotLength/noiseVar;
        Pint = Ue_trans*[1:E]'/slotLength/noiseVar;
        Pint2 = Ue_trans*([1:E].^2)'/slotLength/noiseVar;
        w(2:end) = mean(1-err(P,Pint,Pint2),1);
    end
else
    % transmit power profile of the other devices
    Ue_trans = mnrnd(U-1,pE_transmit,nMC);
    Ue_trans(:,1) = [];
    
    % probability of being in a singleton slot
    p_singleton = sum(sum(Ue_trans,2)==0)/nMC;

    % successful delivery probabilities
    for ii = 1:E
        w(ii+1) = (1-err(ii/noiseVar/slotLength))*p_singleton;
    end
end

if metric_select(1) || metric_select(2)
%% Transition probabilities of the Markov chain (refresh status, battery level)
% See Fig. 3 in [1]
T = zeros(E+1,E+1);
for i = 0:E
for j = 0:E
    if i==j 
        if i  < E, T(i+1,j+1) = (1-eta)*(1-rho(i+1)); end
        if i == E, T(i+1,j+1) = (1-rho(i+1)); end
    elseif j == i + 1
        T(i+1,j+1) = eta*(1-rho(i+1)); 
    elseif j== 0
        T(i+1,j+1) = rho(i+1)*(1-w(i+1)); 
    end
end
end

t0 = rho.*w;

%% Moments of the inter-refresh time Y 
meanY = (eye(E+1)-T)\ones(E+1,1);
meanY2 = meanY + 2*(eye(E+1)-T)^(-2)*T*ones(E+1,1);
meanY = meanY(1);
meanY2 = meanY2(1);
end

%% Average AoI
if metric_select(1)
    avgAoI = 1 + meanY2/2/meanY; 
    output = [output; avgAoI];
end

%% AVP
if metric_select(2)
    pmf_Y = zeros(AoI_thres,1);
    tmp = [1 zeros(1,E)];
    for yy = 1:AoI_thres
        pmf_Y(yy) = tmp*t0;
        tmp = tmp*T;
    end
    
    tmp0 = [1 zeros(1,E)];
    AVP = 1 - ((AoI_thres-1)*tmp0*T^(AoI_thres-1)*ones(E+1,1) + ...
            [1:(AoI_thres-1)]*pmf_Y(1:(AoI_thres-1)))/meanY;

    output = [output; AVP];
end

%% Throughput
if metric_select(3)
    throughput = U*pE_transmit*w;
    output = [output; throughput];
end

%%
if debug
    toc
    fprintf('\n');
    fprintf('throughput:  %.4f\n', throughput )
    fprintf('average AoI: %.4f\n', avgAoI)
    fprintf('AVP:         %.4f\n', AVP)

    fprintf('steady-state battery level dist.: [');
    fprintf('%.4f, ', pE_steady(1:end-1));
    fprintf('%.4f]\n', pE_steady(end));

    fprintf('successful delivery probability:  [');
    fprintf('%.4f, ', w(1:end-1));
    fprintf('%.4f]\n', w(end));
    keyboard
end
