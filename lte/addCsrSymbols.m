function [ X ] = addCsrSymbols(X, toneIndex)
% CSR Symbols 
%
% Two CSR symbols per resource block om OFDM symbols {0, 5, 7, 12} of each
% subframe. There are 8 CSR symbols per resource block, two in each of
% these four OFDM symbols.
%
% Input 
% X             -> FFT Grid
% toneIndex     -> is a look-up table for the tones (in the
%                   FFT) corresponding to each entry of the grid

nSymbols = size(X, 2);
nLayers = size(X, 3);

nREs = length(toneIndex);
nRBs = round(nREs / 12);


% Error checking
if (nSymbols ~= 140)
    error('Invalid number of OFDM symbols in the grid');
end

for iSubframe = 1:10 % For each subframe
    % Generate CSR symbols
    csr = CSRgenerator(iSubframe, nLayers, nRBs);
    
    % Add in the grid
    for iSlot = 1:2 % For each slot
        for iCsrSymbol = 1:2 % 2 symbols with CSR in each slot
            if (iCsrSymbol == 1)
                % CSR Symbols are placed in subcarrier {0,6} of the RB
                iCSRs = toneIndex(1:6:nREs);
                
            else
                % CSR Symbols are placed in subcarrier {3,9} of the RB
                iCSRs = toneIndex(4:6:nREs);
            end
            % Absolute Symbol index within the frame
            iSymbol = (iSubframe - 1)*14 + (iSlot - 1)*7 + ...
                (iCsrSymbol - 1)*4 + 1;
            X(iCSRs, iSymbol, :) = csr(:, iSlot, iCsrSymbol, :);
        end
    end
    
end

