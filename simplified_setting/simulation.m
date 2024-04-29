function [avgAoI,AVP,throughput] = simulation(U,E,alpha,eta,ptx,slotLength,...
            R,noiseVar,capture,SIC,AoI_thres,nSlots)
% function [avgAoI,AVP,throughput] = simulation(U,B,alpha,eta,ptx,slotLength,...
%            R,noiseVar,capture,SIC,AoI_thres,nSlots)
%
% Evaluate the average age of information (AoI), age violation probability 
% (AVP), and throughput of a slotted ALOHA system with energy harvesting 
% via a simulation of the protocol.
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
%   nSlots     : number of simulated slots 
% 
% OUTPUTS:
%   avgAoI : average AoI
%   AVP    : age violation probability
%   throughput : average number of packet delivered per slot

addpath('../helpers')

%% Debug mode?
debug = 1;

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

    capture = 1;
    SIC = 1;

    AoI_thres = 10000;  

    nSlots  = 1e5;              
end

%% Function to compute the error probability in AWGN channel
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

%% Simulation
avgAoI = 0;
AVP = 0;

currAoI = 2*ones(U,1);      % current AoI for each device
nDelivered    = zeros(U,1);       % number of delivered updates by each device

lastUpdate_gen_time = zeros(U,1);   % generation time of the lattest update from each device
last_refresh_time = zeros(U,1);     % lattest refresh time of each device

bat = E*ones(U,1);                % current battery level, assume to be full initially

inter_refresh_time = [];    % inter refresh time, Y

for idxSlot = 1:nSlots
    % update generation 
    new_update = (rand(U,1) <= alpha);
    lastUpdate_gen_time(new_update) = idxSlot;

    % transmission
    trans = zeros(U,1);
    trans(new_update) = (bat(new_update) > 0).*...
        (rand(sum(new_update),1) <= ptx(bat(new_update)+1));
    trans = logical(trans);
    
    sumpower = sum(bat(trans));
    activeUser = find(trans);

    % AoI evolution
    if idxSlot > 1
        currAoI = currAoI + 1;
    end
    
    % content resolution
    if capture || (~capture && sum(trans) == 1) 
        if capture
            w_SIC = compute_w_capture(bat(trans),slotLength,noiseVar,1,err,err_precomp);
        end
        iitmp = 1;
        for idxUser = activeUser(:)'
            if capture 
                decoded = rand <= w_SIC(iitmp);
                iitmp = iitmp + 1;
            else
                decoded = rand <= 1-err(bat(idxUser)/(sumpower - bat(idxUser) ...
                    + slotLength*noiseVar)); 
            end
            if decoded
                nDelivered(idxUser) = nDelivered(idxUser) + 1;
                if idxSlot - lastUpdate_gen_time(idxUser) + 1 < currAoI(idxUser) 
                    inter_refresh_time = [inter_refresh_time idxSlot-last_refresh_time(idxUser)];
                    last_refresh_time(idxUser) = idxSlot;
                end
                currAoI(idxUser) = idxSlot - lastUpdate_gen_time(idxUser) + 1;
            end
        end
    end
    
    % collect the current AoI and AVP values, from slot 101 onward
    if idxSlot > 100
        avgAoI = avgAoI + sum(currAoI);
        AVP = AVP + sum(currAoI > AoI_thres);
    end

    % energy harvesting
    new_energy = (rand(U,1) <= eta);
    bat = min(bat + new_energy, E);
    bat(trans) = 0;
end

%% Compute the average AoI, AVP, and throughput
throughput = sum(nDelivered,'all')/nSlots;
avgAoI = avgAoI/U/(nSlots-100);
AVP = AVP/U/(nSlots-100);

%%
if debug
    toc
    fprintf('\n');
    fprintf('throughput:  sim %.4f\n', throughput)
    fprintf('average AoI: sim %.4f\n', avgAoI)
    fprintf('AVP:         sim %.4f\n', AVP)
    keyboard
end
