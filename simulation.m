function [avgAoI,AVP,S_sim,EZ,EY,EY2,pmf_Z,pmf_Y,pmf_AoI] = ...
            simulation(U,E,alpha,eta,ptx,pretx,slotLength,R,noiseVar,capture, ...
            SIC,AoI_thres,nSlots)

% function [avgAoI,AVP,throughput] = simulation(U,B,alpha,eta,ptx,slotLength,...
%            R,noiseVar,capture,SIC,AoI_thres,nSlots)
%
% Evaluate the average age of information (AoI), age violation probability 
% (AVP), and throughput of a slotted ALOHA system with energy harvesting 
% via a simulation of the protocol.
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
%   nSlots     : number of simulated slots 
% 
% OUTPUTS:
%   avgAoI : average AoI
%   AVP    : age violation probability
%   throughput : average number of packet delivered per slot

addpath('./helpers')

%% Debug mode?
debug = 1;

if debug
    tic
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

    nSlots = 1e5;
end
AoI_max = AoI_thres;

%% Function to compute the error probability in AWGN channel
if capture
    disper = @(P,Pint,Pint2) log2(exp(1))^2.*(P.^2.*(1+2.*Pint + Pint.^2 ...
                - Pint2) + 2*P.*(Pint+1).^3)...
            ./(2*(Pint + 1).^2.*(P+Pint+1).^2);
    err = @(P,Pint,Pint2) qfunc(sqrt(slotLength./disper(P,Pint,Pint2)).*(log2(1+P./(Pint+1))/2 - R));
    % precompute the error prob. in case of no interference (to speed up
    % the computation)
    err_precomp = err((1:E)/slotLength/noiseVar,0,0);
else
    disper = @(SNR) log2(exp(1))^2.*SNR.*(SNR+2)./(1+SNR).^2/2;
    err = @(SNR) qfunc(sqrt(slotLength./disper(SNR)).*(log2(1+SNR)/2 - R));
end

%% Simulation
avgAoI = 0;
AVP = 0;

currAoI = 2*ones(U,1);  % current AoI
nDec    = zeros(U,1);       % number of delivered updates
pmf_AoI = zeros(AoI_max,1);

lastUpdate_gen_time = zeros(U,1);   % generation time of the lattest update
last_refresh_time = zeros(U,1);     % last refresh time

refresh_time = [];
inter_refresh_time = [];

bat = E*ones(U,1);                % current battery level, assume to be full initially

for idxSlot = 1:nSlots
    % update generation and transmission
    new_update = (rand(U,1) <= alpha);
    lastUpdate_gen_time(new_update) = idxSlot;
    trans_power = zeros(U,1);

    % transmit new updates
    bat_tx = new_update.*bat;
    bat_retx = (1-new_update).*bat;
    ntx = sum(bat_tx == (1:E),1);
    nretx = sum(bat_retx == (1:E),1);
    for bb = 1:E
        if ntx(bb)
            trans_power(bat_tx == bb) = randi_distr(0:E, ptx(bb+1,:), ntx(bb), 1);
        end
        if nretx(bb)
            trans_power(bat_retx == bb) = randi_distr(0:E, pretx(bb+1,:), nretx(bb), 1);
        end
    end
    trans = logical(trans_power);
    
    sumpower = sum(trans_power);
    activeUser = find(trans);

    % AoI evolution
    if idxSlot > 1
        currAoI = currAoI + 1;
    end
    
    % content resolution
    if capture || (~capture && sum(trans) == 1)
        if capture
            w_SIC = compute_w_capture(trans_power(trans),slotLength,noiseVar,1,err,err_precomp);
        end
        iitmp = 1;
        for idxUser = activeUser(:)'
            if capture 
                decoded = rand <= w_SIC(iitmp);
                iitmp = iitmp + 1;
            else
                decoded = rand <= 1-err(trans_power(idxUser)...
                    /(sumpower - trans_power(idxUser) + slotLength*noiseVar)); 
            end
            if decoded
                nDec(idxUser) = nDec(idxUser) + 1;
                if idxSlot - lastUpdate_gen_time(idxUser) + 1 < currAoI(idxUser) 
                    refresh_time = [refresh_time idxSlot - lastUpdate_gen_time(idxUser) + 1];
                    inter_refresh_time = [inter_refresh_time idxSlot-last_refresh_time(idxUser)];
                    last_refresh_time(idxUser) = idxSlot;
                end
                currAoI(idxUser) = idxSlot - lastUpdate_gen_time(idxUser) + 1;
            end
        end
    end
    
    % collect the current AoI and AVP values
    avgAoI = avgAoI + sum(currAoI);
    AVP = AVP + sum(currAoI > AoI_thres);

    for iii = 1:U
        if currAoI(iii) <= AoI_max && idxSlot > 100
            pmf_AoI(currAoI(iii)) = pmf_AoI(currAoI(iii)) + 1;
        end
    end

    % energy harvesting
    new_energy = (rand(U,1) <= eta);
    bat = min(bat + new_energy - trans_power, E);
end

%% Compute the average AoI, AVP, and throughput
S_sim = sum(nDec,'all')/nSlots;
avgAoI = avgAoI/U/nSlots; 
AVP = AVP/U/nSlots; 

%% Distribution of Y and Z
EZ = mean(refresh_time);

pmf_AoI = pmf_AoI/(nSlots-100)/U;

EY = mean(inter_refresh_time);
EY2 = mean(inter_refresh_time.^2);

setval = 1:AoI_thres; %sort(unique(inter_refresh_time));
pmf_Y = zeros(size(setval));
for ii = 1:length(setval)
    val = setval(ii);
    pmf_Y(ii) = sum(inter_refresh_time == val);
end
pmf_Y = pmf_Y/length(inter_refresh_time);

pmf_Z = zeros(size(setval));
for ii = 1:length(setval)
    val = setval(ii);
    pmf_Z(ii) = sum(refresh_time == val);
end
pmf_Z = pmf_Z/length(refresh_time);
xlim([1 AoI_thres])

%%
if debug
    toc
    close all
    
    fprintf('\n');
    fprintf('throughput:  %.4f\n', S_sim )
    fprintf('average AoI: %.4f\n', avgAoI)
    fprintf('AVP:         %.4f\n', AVP)
    fprintf('average fresh value: %.4f\n', EZ)
    fprintf('average inter-refresh time: %.4f\n', EY)
    fprintf('average inter-refresh time squared: %.4f\n', mean(inter_refresh_time.^2))
    keyboard
end
end