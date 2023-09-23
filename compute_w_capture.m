function w = compute_w_capture(powers,slotLength,noiseVar,SIC,err,err_precomp)
% function w = compute_w_capture(powers,slotLength,noise,err,SIC,err_precomp)
% 
% Evaluate the successful decoding probability for each device when the
% decoding is with capture
%
% ﻿[1] Khac-Hoang Ngo, G. Durisi, A. Graell i Amat, A. Munari, and F. Lazaro, 
% “Age of information in slotted ALOHA with energy harvesting,” in IEEE 
% Global Communications Conference (Globecom), Kuala Lumpur, Malaysia, Dec. 
% 2023. [Online]. Available: https://research.chalmers.se/publication/537484/file/537484_Fulltext.pdf.
%
% INPUTS: 
%   powers: transmit powers of the devices
%   slotLength : number of channel uses in a slot
%   noiseVar   : noise variance
%   SIC        : indicator if successive interference cancellation is
%                considered (1) or not (0)
%   err        : the function to compute the error probability
%   err_precomp : error probabilities in the case of no interference
% 
% OUTPUTS:
%   w : successful decoding probability

% Sort the devices in ascending order of transmit power
[powers,I] = sort(powers);
[~,J] = sort(I);

% Compute the successful decoding probability
w = zeros(size(powers));
for idxE = length(powers):-1:1
    if powers(idxE) > 0
        if SIC
            if idxE == length(powers)
                P = powers(idxE)/slotLength/noiseVar;
                Pint = sum(powers(1:idxE-1))/slotLength/noiseVar;
                Pint2 = sum(powers(1:idxE-1).^2)/slotLength/noiseVar;
                w(idxE) = 1-err(P,Pint,Pint2);
            elseif powers(idxE) == powers(idxE+1)
                w(idxE) = w(idxE+1);
            elseif powers(idxE) < powers(idxE+1) 
                if idxE == 1
                    w(idxE) = (1-err_precomp(powers(idxE)))*prod(w(idxE+1:end));
                else
                    P = powers(idxE)/slotLength/noiseVar;
                    Pint = sum(powers(1:idxE-1))/slotLength/noiseVar;
                    Pint2 = sum(powers(1:idxE-1).^2)/slotLength/noiseVar;
                    w(idxE) = (1-err(P,Pint,Pint2))*prod(w(idxE+1:end));
                end
            end
        else
            P = powers(idxE)/slotLength/noiseVar;
            Pint = (sum(powers) - powers(idxE))/slotLength/noiseVar;
            Pint2 = (sum(powers.^2) - powers(idxE)^2)/slotLength/noiseVar;
            w(idxE) = 1-err(P,Pint,Pint2);
        end
    end
end

% Sort the devices back in the original order
w = w(J);

