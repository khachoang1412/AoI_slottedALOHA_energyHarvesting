function [avgAoI,EY] = exact(U,E,alpha,eta,ptx,pretx,slotLength,R,noiseVar,capture,SIC)
% function [avgAoI,EY] = exact(U,E,alpha,eta,ptx,pretx,slotLength,R,noiseVar,capture,SIC)
%
% Evaluate the exact average age of information (AoI) of a slotted ALOHA system 
% with energy harvesting via a first-order Markov analysis.
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
% 
% OUTPUTS:
%   avgAoI : average AoI

addpath('./helpers')

%% Debug mode?
debug = 1;

if debug
    tic
    U = 10;           % number of devices
    E = 2;            % battery capacity

    alpha = .1/U;     % update generation rate
    eta = 0.05;       % energy harvesting rate 

    
    ptx =   [1 0 0; 
             0 1 0;
             0 0 1];       % probability of transmitting a fresh update
    pretx = [1 0 0; 
             1 0 0;
             1 0 0];       % probability of retransmitting an old update
    
    slotLength = 100;     % length of a slot
    R = .8;               % transmission rate

    noiseVar = db2pow(-20);

    capture = 0;
    SIC = 1;
end

% This code works for E <= 2 only
if E > 2
    error('The exact analysis currently works for battery capacity at most 2!')
end

% probability that a device accesses the channel
rho = alpha*ptx + (1-alpha)*pretx;

%% Transition probabilities between battery levels of a device
pcond_E = zeros(E+1,E+1);
for i = 0:E
for j = 0:E
    pcond_E(i+1,j+1) = phi(i,j,rho,eta,E);
end
end

%% Generate all possible battery profiles of U - 1 devices
user_bat_profiles_0 = integer_partitions(U-1,E+1);
user_bat_profiles = [];
for ii = 1:size(user_bat_profiles_0,1)
    user_bat_profiles = [user_bat_profiles; uniqueperms(user_bat_profiles_0(ii,:))];
end
n_bat_profiles = size(user_bat_profiles,1);

%% Transition probabilities between battery profiles
Q_bat_profiles = zeros(n_bat_profiles,n_bat_profiles);
for ii1 = 1:n_bat_profiles
for ii2 = 1:n_bat_profiles
    Q_bat_profiles(ii1,ii2) = bat_transitions(user_bat_profiles(ii1,:),user_bat_profiles(ii2,:),pcond_E);
end
end

%% Function to compute error probability
if capture
    disper = @(P,Pint,Pint2) log2(exp(1))^2.*(P.^2.*(1+2.*Pint + Pint.^2 ...
                - Pint2) + 2*P.*(Pint+1).^3)...
            ./(2*(Pint + 1).^2.*(P+Pint+1).^2);
    err = @(P,Pint,Pint2) qfunc(sqrt(slotLength./disper(P,Pint,Pint2)).*(log2(1+P./(Pint+1))/2 - R));
    err_precomp = err((1:E)/slotLength/noiseVar,0,0);
else
    disper = @(SNR) log2(exp(1))^2.*SNR.*(SNR+2)./(1+SNR).^2/2;
    err = @(SNR) qfunc(sqrt(slotLength./disper(SNR)).*(log2(1+SNR)/2 - R));
    err_precomp = err((1:E)/slotLength/noiseVar);
end

%% Successful decoding probability given the transmit power and battery profile of other devices 
w = zeros(E+1,n_bat_profiles);
for ii1 = 1:E
    for ii2 = 1:n_bat_profiles
        w(ii1+1,ii2) = compute_w(ii1,user_bat_profiles(ii2,:),...
                    noiseVar,slotLength,rho,err,capture,SIC,err_precomp);
    end
end

%% Analyse the Markov chain {successfull delivery of the user of interest, 
%   battery level of the user of interest, battery profile of all other users}

% Store all states
nState = 5*(E+1)*n_bat_profiles;
States = [];
for ss = [0.1 0.2 0.3 1 2]    % [0.1 0.2 .3 1 2] corresponds to [S F_new F_old R_new R_old] in the paper
for ee = 0:E
    States = [States; repmat([ss ee],n_bat_profiles,1) (1:n_bat_profiles)'];
end
end

% Transition probabilities between the states
Q = zeros(nState,nState);
for ii1 = 1:nState
for ii2 = 1:nState
    if States(ii1,1) == .1 %|| States(ii1,1) == 1 || States(ii1,1) == 2
        if States(ii2,1) == 1
            Q(ii1,ii2) = Q_bat_profiles(States(ii1,3),States(ii2,3))*alpha*...
                phi_w(States(ii1,2),States(ii2,2),ptx,eta,E,w,States(ii1,3));
        elseif States(ii2,1) == .1
            Q(ii1,ii2) = Q_bat_profiles(States(ii1,3),States(ii2,3))*(1-alpha)*...
                phi(States(ii1,2),States(ii2,2),pretx,eta,E);
        elseif States(ii2,1) == .2
            Q(ii1,ii2) = Q_bat_profiles(States(ii1,3),States(ii2,3))*alpha*...
                phi_w(States(ii1,2),States(ii2,2),ptx,eta,E,1-w,States(ii1,3));
        end
    elseif States(ii1,1) == .2 %|| States(ii1,1) == .3
        if States(ii2,1) == 1 
            Q(ii1,ii2) = Q_bat_profiles(States(ii1,3),States(ii2,3))*alpha*...
                phi_w(States(ii1,2),States(ii2,2),ptx,eta,E,w,States(ii1,3));
        elseif States(ii2,1) == 2 
            Q(ii1,ii2) = Q_bat_profiles(States(ii1,3),States(ii2,3))*(1-alpha)*...
                phi_w(States(ii1,2),States(ii2,2),pretx,eta,E,w,States(ii1,3));
        elseif States(ii2,1) == .2 
            Q(ii1,ii2) = Q_bat_profiles(States(ii1,3),States(ii2,3))*alpha*...
                phi_w(States(ii1,2),States(ii2,2),ptx,eta,E,1-w,States(ii1,3));
        elseif States(ii2,1) == .3
            Q(ii1,ii2) = Q_bat_profiles(States(ii1,3),States(ii2,3))*(1-alpha)*...
                phi_w(States(ii1,2),States(ii2,2),pretx,eta,E,1-w,States(ii1,3));
        end
    end
end
end

blocksize = nState/5;

% Q = Q./(sum(Q,2));

Q(2*blocksize+1:3*blocksize,:) = Q(blocksize+1:2*blocksize,:);
Q(3*blocksize+1:4*blocksize,:) = Q(1:blocksize,:);
Q(4*blocksize+1:5*blocksize,:) = Q(1:blocksize,:);

%% Compute the mean of Y
Qss = Q(1:blocksize,1:blocksize);
Qsf = Q(1:blocksize,blocksize+1:3*blocksize);
rs = sum(Q(1:blocksize,3*blocksize+1:5*blocksize),2);
Qff = Q(blocksize+1:3*blocksize,blocksize+1:3*blocksize);
rf = sum(Q(blocksize+1:3*blocksize,3*blocksize+1:nState),2);

ef = (eye(2*blocksize) - Qff)\(1+rf);
es = (eye(blocksize) - Qss)\(1+rs+Qsf*ef);
er = ones(1,2*blocksize);

weights = sum(Q(3*blocksize+1:end,:));
weights = weights'/sum(weights);

EY = [es' ef' er]*weights;

%% Compute the mean of Y^2
ef2 = (eye(2*blocksize) - Qff)\(1 + 3*rf + 2*Qff*ef);
es2 = (eye(blocksize) - Qss)\(1 + 3*rs + 2*Qsf*ef + 2*Qss*es + Qsf*ef2);
er2 = ones(1,2*blocksize);

EY2 = [es2' ef2' er2]*weights;


%% Steady state distribution
Qhat = Q';

Qhat(1,:) = 1;
for ii = 2:nState
    Qhat(ii,ii) = Qhat(ii,ii) - 1;
end

Qhat_inv = pinv(Qhat);
p_steady = Qhat_inv(:,1); 
p_steady = max(p_steady,0);
p_steady = p_steady/sum(p_steady);

%% Find backward transition probabilities
Q_backward = ((Q./p_steady').*p_steady)';
Q_backward(isnan(Q_backward)) = 0;
Q_backward(Q_backward == inf) = 0;
Qff_backward = Q_backward(2*blocksize+1:3*blocksize,2*blocksize+1:3*blocksize);

%% Compute the mean of Z
rfZ = sum(Q_backward(2*blocksize+1:3*blocksize,[1:2*blocksize 3*blocksize+1:nState]),2);
efZ = (eye(blocksize) - Qff_backward)\(1 + rfZ);

weights = Q_backward(4*blocksize+1:5*blocksize,...
    [blocksize+1:2*blocksize 2*blocksize+1:3*blocksize])';
% weights = weights'/sum(weights;

e2Z = [ones(length(efZ),1); efZ]'*weights;

r1Z = sum(p_steady(3*blocksize+1:4*blocksize));
weight2Z = p_steady(4*blocksize+1:5*blocksize);
r12Z = r1Z + sum(sum(weights)*weight2Z);
EZ = (r1Z + e2Z*weight2Z)/r12Z;

%% Average AoI
avgAoI = EZ + EY2/2/EY;

%%
if debug
    toc
    close all
    
    fprintf('\n');
    fprintf('average AoI: %.4f\n', avgAoI)
    fprintf('average fresh AoI value: %.4f\n', EZ)
    fprintf('average inter-refresh time: %.4f\n', EY)
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
function out = phi_w(i,j,pp,eta,B,w,idxBatProfile)
    out = 0;
    if i >= j
        out = out + pp(i+1,i-j+1)*(1-eta*((i~=B) || (j~=B)))*w(i-j+1,idxBatProfile);
    end
    if i >= j-1 && i-j+1 <= B
        out = out + pp(i+1,i-j+2)*eta*w(i-j+2,idxBatProfile);
    end
end

%% Compute transition probabilities from battery profile Upre to Unext
% This function works for Psi <= 2 only
function [ptrans,uu,scale] = bat_transitions(Upre,Unext,pcond_E)
    uu = [];
    scale = [];
    if size(pcond_E,1) == 2
        for u10 = 0:Upre(2)
            u11 = Upre(2) - u10;
            u00 = Unext(1) - u10;
            u01 = Upre(1) - u00;
            if sum([u00 u10 u01 u11] < 0) == 0
                uu = [uu; u00 u10 u01 u11];
                scale = [scale; nchoosek(Upre(1),u00)*nchoosek(Upre(2),u10)];
            end
        end
    elseif size(pcond_E,1) == 3
        u02 = 0;
        for u00 = 0:floor(min(Upre(1),Unext(1)))
            u01 = Upre(1) - u00;
            for u10 = 0:floor(min(Upre(2),Unext(1)-u00))
                u20 = Unext(1) - u00 - u10;
                for u11 = 0:floor(min([Upre(2)-u10 Unext(2)-u01]))
                    u21 = Unext(2) - u01 - u11;
                    u22 = Upre(3) - u20 - u21;
                    u12 = Unext(3) - u22;
                    if sum([u00 u10 u20 u01 u11 u21 u02 u12 u22] < 0) == 0
                        uu = [uu; u00 u10 u20 u01 u11 u21 u02 u12 u22];
                        scale = [scale; nchoosek(Upre(1),u00)*...
                            nchoosek(Upre(2),u10)*nchoosek(Upre(2)-u10,u11)*...
                                nchoosek(Upre(3),u20)*nchoosek(Upre(3)-u20,u21)];
                    end
                end
            end
        end
    end
    if sum(uu(:) < 0)
        keyboard
    end

    ptrans = 0;
    if ~isempty(uu)
        ptrans = sum(prod(pcond_E(:)'.^uu,2).*scale);
    end
end

%% Compute the successful decoding probability
function w0 = compute_w(E0,Ue,noise,slotLength,rho,err,capture,SIC,err_precomp)
    E = length(Ue)-1;
    w0 = 0;
    if E0 > 0
        if capture 
            nMC = 1e4;
            Ue_trans = zeros(nMC,E,E+1);
            for idxMC = 1:nMC
                for idxB = 1:E
                    Ue_trans(idxMC,idxB,:) = mnrnd(Ue(idxB+1),rho(idxB+1,:),1); 
                end
            end
            Ue_trans = squeeze(sum(Ue_trans,2)); Ue_trans(:,1) = [];
            
            if SIC
                w0_tmp = 0;
                for idxMC = 1:nMC
                        powers = [E0 flip(repelem(1:E,Ue_trans(idxMC,:)))];
                        w_tmp = compute_w_capture(powers,slotLength,noise,1,err,err_precomp);
                            w0_tmp = w0_tmp + w_tmp(1); 
                end
                w0 = w0_tmp/nMC;
            else
                P = E0/slotLength/noise;
                Pint = Ue_trans*[1:E]'/slotLength/noise;
                Pint2 = Ue_trans*([1:E].^2)'/slotLength/noise;
                w0 = mean(1-err(P,Pint,Pint2),1);
            end
        else
            w0 = (1-err(E0/noise/slotLength))*prod(rho(:,1).^Ue(:));
        end
    end
end