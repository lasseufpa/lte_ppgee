function [ FEQ ] = feqFromCsrSymbols(X, Y)
% Estimate FEQ from CSR Symbols
%
% This function extracts the CSRs from the Rx grid and compares to the CSRs
% that should be obtained (as in the Tx grid) to the obtain the FEQ and
% sparse points in the grid. It, then, interpolates the FEQ to obtain the
% equalizer gain in all Resource Elements. In the end, one FEQ matrix is
% output with number of rows equivalent to the number of resource elements
% and number of columns equivalent to the number of subframes in a frame.
% The latter is explained by the fact that one distinct FEQ is applied in
% each subframe.
%
% Note two CSR symbols per resource block om OFDM symbols {0, 4, 7, 11} (or
% {1, 5, 8, 12} in MATLAB) of each subframe. There are 8 CSR symbols per
% resource block, two in each of these four OFDM symbols.
%
% Input
% X             -> Grid of Tx Symbols
% Y             -> Grid of Rx Symbols
%
% Output
% FEQ           -> Matrix with FEQ taps (dimensions of "nREs x 10")


% Infer Parameters
nREs     = size(Y, 1);
nSymbols = size(Y, 2);
nLayers  = size(Y, 3);
nRBs    = round(nREs / 12);

% Preallocate
H = zeros(nRBs*4, 10);
H_grid = zeros(nREs, 10);

% Error checking
if (nSymbols ~= 140)
    error('Invalid number of OFDM symbols in the grid');
end


%% Estimate channel response

% If we look anternately over the symbols that carry CSRs, they alternated
% in the RE indexes where the pilots are placed. Thus, we form two sets of
% indexes and use those to extract exactly the CSRs.
iReCSRsEven = 1:6:nREs;
iReCSRsOdd  = 4:6:nREs;
iReCSRs     = union(iReCSRsEven, iReCSRsOdd);

for iSubframe = 1:10 % For each subframe

    % Absolute Symbol index within the frame 
    iEvenSymbol = (iSubframe - 1)*14 + [1 8];
    iOddSymbol  = (iSubframe - 1)*14 + [5 12];

    % Channel tap estimation:
    H(1:2:end, iSubframe) = sum(Y(iReCSRsEven, iEvenSymbol) ...
        ./ X(iReCSRsEven, iEvenSymbol), 2) / 2;
    H(2:2:end, iSubframe) = sum(Y(iReCSRsOdd, iOddSymbol) ...
        ./ X(iReCSRsOdd, iOddSymbol), 2) / 2;

end

%% Interpolate
% There are only 4 Resource Elements per Resource Block with CSRs. Thus,
% the number of channel taps that were estimated for each subframe and RB
% was 4. We need now to fill all the 12 REs of each RB with estimated taps.
% This is done by interpolation.

H_grid(iReCSRs, :) = H;

L_mov = 4;
h_mov = (1/L_mov)*ones(L_mov, 1);

%% Obtain FEQ

FEQ = 1./H_grid;
end
