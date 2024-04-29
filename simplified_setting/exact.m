function avgAoI = exact(U,E,alpha,eta,ptx,slotLength,R,noiseVar,capture,SIC)
% function avgAoI = exact(U,B,alpha,eta,ptx,slotLength,R,noiseVar,capture,SIC)
%
% Evaluate the exact average age of information (AoI) of a slotted ALOHA system 
% with energy harvesting via a first-order Markov analysis.
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
% 
% OUTPUTS:
%   avgAoI : average AoI

addpath('../helpers')

%% Debug mode?
debug = 0;

if debug
    tic
    U = 10;         
    E = 2;  

    alpha = 1/U;      
    eta = 0.05;                   
    ptx = [0 1 1]'; 
    
    slotLength = 100;      
    R = .8;             
    noiseVar = db2pow(-20);

    capture = 1;
    SIC = 1;
end

% This code works for E <= 2 only
if E > 2
    error('The exact analysis currently works for battery capacity at most 2!')
end

% probability that a device accesses the channel (i.e., generates a new 
% update and transmits it) 
rho = alpha*ptx;

%% Transition probabilities between battery levels of a device
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

%% Generate all possible battery profiles of U - 1 devices
bat_profiles_0 = integer_partitions(U-1,E+1);
bat_profiles = [];
for ii = 1:size(bat_profiles_0,1)
    bat_profiles = [bat_profiles; uniqueperms(bat_profiles_0(ii,:))];
end
n_bat_profiles = size(bat_profiles,1);

%% Transition probabilities between battery profiles
Q_bat_profiles = zeros(n_bat_profiles,n_bat_profiles);
for ii1 = 1:n_bat_profiles
    for ii2 = 1:n_bat_profiles
        Q_bat_profiles(ii1,ii2) = bat_transitions(bat_profiles(ii1,:),...
            bat_profiles(ii2,:),pE_transition);
    end
end

%% Function to compute error probability
if capture
    dispersion = @(P,Pint,Pint2) log2(exp(1))^2.*(P.^2.*(1+2.*Pint + Pint.^2 ...
                - Pint2) + 2*P.*(Pint+1).^3)...
            ./(2*(Pint + 1).^2.*(P+Pint+1).^2);
    err = @(P,Pint,Pint2) qfunc(sqrt(slotLength./dispersion(P,Pint,Pint2))...
        .*(log2(1+P./(Pint+1))/2 - R));
    err_precomp = err((1:E)/slotLength/noiseVar,0,0);
else
    dispersion = @(SNR) log2(exp(1))^2.*SNR.*(SNR+2)./(1+SNR).^2/2;
    err = @(SNR) qfunc(sqrt(slotLength./dispersion(SNR)).*(log2(1+SNR)/2 - R));
    err_precomp = err((1:E)/slotLength/noiseVar);
end

%% Successful decoding probability given the transmit power and battery profile of other devices 
w = zeros(E+1,n_bat_profiles);
for ii1 = 1:E
    for ii2 = 1:n_bat_profiles
        w(ii1+1,ii2) = compute_w(ii1,bat_profiles(ii2,:),...
                    noiseVar,slotLength,ptx,alpha,err,capture,SIC,err_precomp);
    end
end

%% Analyse the Markov chain {successfull delivery of the user of interest, 
%   battery level of the user of interest, battery profile of all other users}

% Store all states
nState = 2*(E+1)*n_bat_profiles;
States = [];
for ss = 0:1
for ee = 0:E
    States = [States; repmat([ss ee],n_bat_profiles,1) (1:n_bat_profiles)'];
end
end

% Transition probabilities between the states
Q = zeros(nState,nState);
for ii1 = 1:nState
for ii2 = 1:nState
    if States(ii1,1) == 1 && States(ii2,1) == 0 
        Q(ii1,ii2) = Q_bat_profiles(States(ii1,3),States(ii2,3))*...
            ((1-eta)*(States(ii2,2) == 0) + eta*(States(ii2,2) == 1));
    elseif States(ii1,1) == 0 && States(ii2,1) == 0
        Q(ii1,ii2) = Q_bat_profiles(States(ii1,3),States(ii2,3))*...
            ( (1-alpha*ptx(States(ii1,2)+1)*(States(ii1,2) > 0))*...
                ((1-eta)*(States(ii2,2) == States(ii1,2) && States(ii1,2) < E) +...
                    eta*(States(ii2,2) == States(ii1,2)+1) + ...
                    (States(ii2,2) == States(ii1,2) && States(ii1,2) == E)) +...
            alpha*ptx(States(ii1,2)+1)*(States(ii1,2) > 0)*(States(ii2,2) == 0) *...
                (1-w(States(ii1,2)+1,States(ii1,3))));
    elseif States(ii1,1) == 0 && States(ii2,1) == 1 && States(ii1,2) > 0 && States(ii2,2) == 0
        Q(ii1,ii2) = Q_bat_profiles(States(ii1,3),States(ii2,3))*...
            (alpha*ptx(States(ii1,2)+1)*w(States(ii1,2)+1,States(ii1,3)));
    end 
end
end
Q = Q./(sum(Q,2));

%% Compute the mean of the inter refresh time, Y
b = 1+sum(Q(1:nState/2,nState/2+1:end),2);
A = eye(nState/2) - Q(1:nState/2,1:nState/2);
meanY0 = A\b;

weights = sum(Q(nState/2:nState/2+n_bat_profiles,1:2*n_bat_profiles)); 
weights = weights'/sum(weights);

meanY = meanY0(1:2*n_bat_profiles)'*weights;

%% Compute the mean of Y^2
b = b + 2*(meanY0 - 1);
meanY1 = A\b;
meanY2 = meanY1(1:2*n_bat_profiles)'*weights;

%% Average AoI
avgAoI = 1 + meanY2/2/meanY;

%%
if debug
    toc
    fprintf('average AoI: sim %.4f\n', avgAoI)
    keyboard
end

end

%% Compute transition probabilities from battery profile Upre to Unext
% See Eq.(8) in [1]
% This function works for Psi <= 2 only
function [ptrans,nn,scale] = bat_transitions(Upre,Unext,pE_transition)
    nn = [];
    scale = [];
    if size(pE_transition,1) == 2
        for n10 = 0:Upre(2)
            n11 = Upre(2) - n10;
            n01 = Upre(1) - Unext(1) + n10;
            n00 = Upre(1) - n01;
            if sum([n00 n10 n01 n11] < 0) == 0
                nn = [nn; n00 n10 n01 n11];
                scale = [scale; nchoosek(Upre(1),n00)*nchoosek(Upre(2),n10)];
            end
        end
    elseif size(pE_transition,1) == 3
        for n20 = max(ceil(Upre(3) - Unext(3)),0):floor(min(Upre(3),Unext(1)))
            n12 = Unext(3) - Upre(3) + n20;
            for n10 = max(ceil(Upre(2)-Unext(2)-n12),0):floor(min([Upre(2)-n12; Unext(1); -Unext(2)+Upre(2)-n12+Upre(1)]))
                n01 = Unext(2) - Upre(2) + n12 + n10;
                n00 = Upre(1) - n01;
                n11 = Upre(2) - n10 - n12;
                n22 = Upre(3) - n20;
                    nn = [nn; n00 n10 n20 n01 n11 0 0 n12 n22];
                    scale = [scale; nchoosek(Upre(1),n00)*nchoosek(Upre(2),n10)*nchoosek(Upre(2)-n10,n11)*...
                                nchoosek(Upre(3),n20)];
            end
        end
    end

    ptrans = 0;
    if ~isempty(nn)
        ptrans = sum(prod(pE_transition(:)'.^nn,2).*scale);
    end
end

%% Compute the successful decoding probability
function w0 = compute_w(E0,Ue,noiseVar,slotLength,ptx,alpha,err,capture,SIC,err_precomp)
    E = length(Ue)-1;
    w0 = 0;
    if E0 > 0
        if capture 
            nMC = 1e4;
            Ue_trans = zeros(nMC,E);
            for idxMC = 1:nMC
                for idxE = 1:E
                    Ue_trans(idxMC,idxE) = sum(rand(Ue(idxE+1),1) <= alpha*ptx(idxE+1));
                end
            end
            
            if SIC
                w0_tmp = 0;
                for idxMC = 1:nMC
                    powers = [E0 flip(repelem(1:E,Ue_trans(idxMC,:)))];
                    w_tmp = compute_w_capture(powers,slotLength,noiseVar,1,err,err_precomp);
                        w0_tmp = w0_tmp + w_tmp(1); 
                end
                w0 = w0_tmp/nMC; 
            else
                P = E0/slotLength/noiseVar;
                Pint = Ue_trans*[1:E]'/slotLength/noiseVar;
                Pint2 = Ue_trans*([1:E].^2)'/slotLength/noiseVar;
                w0 = mean(1-err(P,Pint,Pint2),1);
            end
        else
            w0 = (1-err(E0/noiseVar/slotLength))*...
                    prod((1-alpha*ptx(:)).^Ue(:));
        end
    end
end