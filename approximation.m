function [output,w,pE_steady] = approximation(U,E,alpha,eta,ptx,pretx, ...
        slotLength,R,noiseVar,capture,SIC,AoI_thres,metric_sel,nMC)
% function [output,w,pE_steady] = approximation(U,alpha,eta,E,...
%    slotLength,R,ptx,pretx,noise,capture,SIC,AoI_thres,metric_sel,nMC)
%
% Evaluate approximations of the average age of information (AoI), age ...
% violation probability (AVP), and throughput of a slotted ALOHA system 
% with energy harvesting.
%
% ﻿Reference:
% 
% ﻿Khac-Hoang Ngo, G. Durisi, A. Munari, and F. Lazaro, and A. Graell i Amat,
% "Timely status updates in slotted ALOHA networks with energy harvesting," 
% submitted to IEEE Transactions on Communications, Apr. 2024.
%
% This code implements a more general version of the protocol considered in
% this reference. Here, we assume that a device can also retransmit
% nonfresh updates.
%
% INPUTS: 
%   U     : number of devices
%   E     : battery capacity of a device
%   alpha : probability that a device has a new update in a slot
%   eta   : probability that a device harvests an energy unit in a slot
%   ptx   : (E+1)x(E+1) matrix whose (i,j)th entry is the 
%           probability that a device transmits with j energy units if it 
%           has a new update and its battery level is i \in [0,1,...,E].
%   pretx : (E+1)x(E+1) matrix whose (i,j)th entry is the 
%           probability that a device retransmits the latest update with j 
%           energy units if its battery level is i \in [0,1,...,E].
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

addpath('./helpers')

%% Debug mode?
debug = 1;

if debug
    U = 1000;           % number of devices
    E = 8;              % battery capacity

    alpha = 2/U;        % update generation rate 
    eta = .005;         % energy harvesting rate
    
    ptx = eye(E+1);
    pretx = zeros(E+1,E);
    pretx = [1-sum(pretx,2) pretx];
    
    slotLength = 100;   % length of a slot
    R = .8;             % transmission rate
    
    AoI_thres = 10000;  % age threshold
    nSlots  = 1e5;      % number of slots

    noiseVar = db2pow(-20);
    
    capture = 1; 
    SIC = 1;

    metric_sel = [1 1 1];
    nMC = 1e5;      % number of Monte-Carlo simulation for the power profile of all devices
end

output = [];

% probability that a device accesses the channel
rho = alpha*ptx + (1-alpha)*pretx;

%% Battery evolution
% --- Transition probabilities between battery levels, pE_transition
pcond_E = zeros(E+1,E+1);
for i = 0:E
for j = 0:E
    pcond_E(i+1,j+1) = phi(i,j,rho,eta,E);
end
end
Qhat = pcond_E';

% --- Steady-state PMF of the battery level, pE_steady
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
    disper = @(P,Pint,Pint2) log2(exp(1))^2.*(P.^2.*(1+2.*Pint + Pint.^2 ...
                - Pint2) + 2*P.*(Pint+1).^3)...
            ./(2*(Pint + 1).^2.*(P+Pint+1).^2);
    err = @(P,Pint,Pint2) qfunc(sqrt(slotLength./disper(P,Pint,Pint2)).*(log2(1+P./(Pint+1))/2 - R));
    % precompute the error prob. in case of no interference (to speed up
    % the computation)
    err_precomp = err((1:E)/slotLength/noiseVar,0,0);
elseif noiseVar > 0
    disper = @(SNR) log2(exp(1))^2.*SNR.*(SNR+2)./(1+SNR).^2/2;
    err = @(SNR) qfunc(sqrt(slotLength./disper(SNR)).*(log2(1+SNR)/2 - R));
end


% --- Compute the successful delivery probability corresponding to each 
% battery level, which involves an expectation over a multinomial distribution

% initialize the successfull delivery probability 
w = zeros(E+1,1);

% now, compute the successfull delivery probability, w
if capture
    if SIC % decoding with capture and SIC
        % generate realizations of the power profile of the active devices 
        Ue_trans = mnrnd(U,pE_steady'*rho,nMC);
        Ue_trans(:,1) = [];
        
        % compute the success probability for each realization
        w_SIC_tmp = zeros(E,1);
        for idxMC = 1:nMC
            if sum(Ue_trans(idxMC,:))
                % transmit powers of each active devices
                powers = flip(repelem(1:E,Ue_trans(idxMC,:)));

                % compute the success probability of each active device
                if isscalar(powers)
                    w_tmp = 1-err_precomp(powers);
                else
                    w_tmp = compute_w_capture(powers,slotLength,noiseVar,1,err,err_precomp);
                end
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
        Ue_trans = mnrnd(U-1,pE_steady'*rho,nMC);
        Ue_trans(:,1) = [];

        % successful delivery probabilities
        P = [1:E]/slotLength/noiseVar;
        Pint = Ue_trans*[1:E]'/slotLength/noiseVar;
        Pint2 = Ue_trans*([1:E].^2)'/slotLength/noiseVar;

        w(2:end) = mean(1-err(P,Pint,Pint2),1);
    end

elseif noiseVar > 0
    % option 1: exact, slow
%     Ue_trans = mnrnd(U-1,pE_transmit,nMC);
%     Ue_trans(:,1) = [];
%     p_singleton = sum(sum(Ue_trans,2)==0)/nMC;
    
    % option 2: very accurate approx., fast
    p_singleton = (pE_steady'*rho(:,1))^(U-1);  
    
    % successful delivery probabilities
    for energy = 1:E
        w(energy+1) = (1-err(energy/noiseVar/slotLength))*p_singleton;
    end

elseif noiseVar == 0  % collision channel
    p_singleton = (pE_steady'*rho(:,1))^(U-1);
    w(2:E+1) = ones(E,1)*p_singleton;
end


%% Transition probabilities of the Markov chain 
if metric_sel(1) + metric_sel(2) > 0
    % ----- Y -----

    % Markov chain X_B
    nB = E+1;
    P = zeros(3*nB,3*nB);     % in the order: R, S, F
    for i = 1:3*nB
    for j = 1:3*nB
        if i <= 2*nB   % starting from Rout or S
            if j >= nB + 1 && j <= 2*nB   % the destination is S 
                P(i,j) = (1-alpha)*phi(rem(i-1,nB),rem(j-1,nB),pretx,eta,E);
            elseif j >= 2*nB + 1 && j <= 3*nB % the destination is F
                P(i,j) = alpha*phi_w(rem(i-1,nB),rem(j-1,nB),ptx,eta,E,1-w);
            elseif j <= nB  % the destination is R
                P(i,j) = alpha*phi_w(rem(i-1,nB),rem(j-1,nB),ptx,eta,E,w);
            end
        elseif i <= 3*nB    % starting from F
            if j >= 2*nB + 1 && j <= 3*nB % the destination is F
                P(i,j) = phi_w(rem(i-1,nB),rem(j-1,nB),rho,eta,E,1-w);
            elseif j <= nB  % the destination is R
                P(i,j) = phi_w(rem(i-1,nB),rem(j-1,nB),rho,eta,E,w);
            end
        end
    end
    end
    
    % steady state distribution
    Qhat = P';
    Qhat(1,:) = 1;
    for b = 2:height(P)
        Qhat(b,b) = Qhat(b,b) - 1;
    end
    Qhat_inv = pinv(Qhat);
    psteady = max(Qhat_inv(:,1),0); 
%     psteady = psteady/sum(psteady);
    pR = psteady(1:nB); pR = pR(:)/sum(pR);

    % split R into Rout and Rin
    nB = E+1;
    P = zeros(4*nB,4*nB);     % in the order: Rout, S, F, Rin
    for i = 1:4*nB
    for j = 1:4*nB
        if i <= 2*nB   % starting from Rout or S
            if j >= nB + 1 && j <= 2*nB   % the destination is S 
                P(i,j) = (1-alpha)*phi(rem(i-1,nB),rem(j-1,nB),pretx,eta,E);
            elseif j >= 2*nB + 1 && j <= 3*nB % the destination is F
                P(i,j) = alpha*phi_w(rem(i-1,nB),rem(j-1,nB),ptx,eta,E,1-w);
            elseif j >= 3*nB + 1  % the destination is Rin
                P(i,j) = alpha*phi_w(rem(i-1,nB),rem(j-1,nB),ptx,eta,E,w);
            end
        elseif i <= 3*nB    % starting from F
            if j >= 2*nB + 1 && j <= 3*nB % the destination is F
                P(i,j) = phi_w(rem(i-1,nB),rem(j-1,nB),rho,eta,E,1-w);
            elseif j >= 3*nB + 1  % the destination is Rin
                P(i,j) = phi_w(rem(i-1,nB),rem(j-1,nB),rho,eta,E,w);
            end
        elseif i==j && i >= 3*nB+1
            P(i,j) = 1;
        end
    end
    end

    % group all the Rin
    TY = P(1:3*nB,1:3*nB);       % between transient states
    t0Y = sum(P(1:3*nB,3*nB+1:end),2);  % from transient states to absorbing states
    
    % moments of Y
    meanY = (eye(3*nB)-TY)\ones(3*nB,1);
    meanY = pR'*meanY(1:nB);
    
    % ----- Z -----
    if sum(pretx(:,2:end),'all')
        Ptmp = zeros(4*nB,4*nB);
        for i = 1:4*nB
        for j = 1:4*nB
            if i <= 2*nB   % starting from S, take alpha = 1
                if j >= 2*nB + 1 && j <= 3*nB % the destination is F
                    Ptmp(i,j) = phi_w(rem(i-1,nB),rem(j-1,nB),ptx,eta,E,1-w);
                elseif j >= 3*nB + 1  % the destination is Rin
                    Ptmp(i,j) = phi_w(rem(i-1,nB),rem(j-1,nB),ptx,eta,E,w);
                end
            elseif i <= 3*nB    % starting from F, take alpha = 0
                if j >= 2*nB + 1 && j <= 3*nB % the destination is F
                    Ptmp(i,j) = phi_w(rem(i-1,nB),rem(j-1,nB),pretx,eta,E,1-w);
                elseif j >= 3*nB + 1  % the destination is Rin
                    Ptmp(i,j) = phi_w(rem(i-1,nB),rem(j-1,nB),pretx,eta,E,w);
                end
            end
        end
        end
        
        TZ = cell(nB,1); % transition between transient states when the 
                        %  device generates an update at battery level b
        t0Z = cell(nB,1); % transition to absorbing states

        normalizeZ = 0;
        for b = 0:E
            TZ{b+1} = zeros(nB+1,nB+1);
            TZ{b+1}(:,2:end) = Ptmp([nB+b+1 2*nB+1:3*nB],2*nB+1:3*nB);

            t0Z{b+1} = sum(Ptmp([nB+b+1 2*nB+1:3*nB],3*nB+1:4*nB),2);

            normalizeZ = normalizeZ + alpha*pE_steady(b+1)*[1 zeros(1,nB)]*((eye(nB+1)-(1-alpha)*TZ{b+1})\t0Z{b+1});
        end
    end
end

%% Average AoI
if metric_sel(1)
    % mean of Z
    if sum(pretx(:,2:end),'all')
        meanZ_num = 0;
        for b = 0:E
            meanZ_num = meanZ_num + alpha*pE_steady(b+1)*[1 zeros(1,nB)]*((eye(nB+1)-(1-alpha)*TZ{b+1})^(-2)*t0Z{b+1});
        end
        meanZ = meanZ_num/normalizeZ;
    else
        meanZ = 1;
    end

    % mean of Y^2
    meanY2 = (eye(3*nB)-TY)^(-2)*ones(3*nB,1);
    meanY2 = 2*pR'*meanY2(1:nB) - meanY;

    % average AoI
    avgAoI = meanZ + meanY2/2/meanY; 

    output = [output; metric_sel(1)*avgAoI];
end

%% Age violation probability

if metric_sel(2)
    % PMF of Y
    pmf_Y = zeros(AoI_thres,1);
    tmp = [pR' zeros(1,2*nB)];
    for yy = 1:length(pmf_Y)
        pmf_Y(yy) = tmp*t0Y;
        tmp = tmp*TY;
    end
    cdf_Y = cumsum(pmf_Y);
    
    
    % PMF of Z
    if sum(pretx(:,2:end),'all')
        pmf_Zb = zeros(AoI_thres,E+1);
        for b = 0:E
            tmp = [1 zeros(1,E+1)];
            for zz = 1:length(pmf_Zb)
                pmf_Zb(zz,b+1) = tmp*t0Z{b+1};
                tmp = tmp*TZ{b+1};
            end
        end
        pmf_Z = pmf_Zb*pE_steady;
        pmf_Z = pmf_Z*alpha.*(1-alpha).^((0:length(pmf_Zb)-1)')/normalizeZ;
    else 
        pmf_Z = [1;zeros(AoI_thres-1,1)];
    end
    
    % Compute AVP
    if sum(pretx(:,2:end),'all')
        AVP = 1 - sum(pmf_Z(1:AoI_thres-1));
        for zz = 1:AoI_thres-1
            AVP = AVP + pmf_Z(zz)*(1 - ((AoI_thres-zz)*(1-cdf_Y(AoI_thres-zz)) + ...
                [1:(AoI_thres-zz)]*pmf_Y(1:(AoI_thres-zz)))/meanY );
        end
    else
        AVP = 1 - ((AoI_thres-1)*(1-cdf_Y(AoI_thres-1)) + ...
                [1:(AoI_thres-1)]*pmf_Y(1:(AoI_thres-1)))/meanY ;
    end

    output = [output; metric_sel(2)*AVP];
end

%% Throughput
if metric_sel(3) 
    S = U*pE_steady'*rho*w;
    output = [output; metric_sel(3)*S];
end


%%
if debug
    toc
    close all
% 
    fprintf('\n');
    if metric_sel(3), fprintf('throughput:  %.4f\n', S); end
    if metric_sel(1), fprintf('average AoI: %.4f\n', avgAoI); end
    if metric_sel(2), fprintf('AVP:         %.4f\n', AVP); end

    fprintf('steady-state battery level dist.: [');
    fprintf('%.4f, ', pE_steady(1:end-1));
    fprintf('%.4f]\n', pE_steady(end));

    fprintf('successful delivery probability:  [');
    fprintf('%.4f, ', w(1:end-1));
    fprintf('%.4f]\n', w(end));
    
    keyboard
end
end

%% Function to compute the probability of transition from battery level i to j
function out = phi(i,j,pp,eta,B)
    out = 0;
    if i >= j
        out = out + pp(i+1,i-j+1)*(1-eta*((i~=B) || (j~=B)));
    end
    if i >= j-1 && i-j+1 <= B
        out = out + pp(i+1,i-j+2)*eta;
    end
end

%% Function to compute the probability of transition from battery level i to j ...
%% and sucessfully/unsucessfully deliver an update
function out = phi_w(i,j,pp,eta,B,w)
    out = 0;
    if i >= j
        out = out + pp(i+1,i-j+1)*(1-eta*((i~=B) || (j~=B)))*w(i-j+1);
    end
    if i >= j-1 && i-j+1 <= B
        out = out + pp(i+1,i-j+2)*eta*w(i-j+2);
    end
end
