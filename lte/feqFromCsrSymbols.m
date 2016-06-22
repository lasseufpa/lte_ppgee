function [ FEQ ] = feqFromCsrSymbols(X, Y)
% CSR Symbols
%
% Two CSR symbols per resource block om OFDM symbols {0, 4, 7, 11} (or {1,
% 5, 8, 12} in MATLAB) of each subframe. There are 8 CSR symbols per
% resource block, two in each of these four OFDM symbols.
%
% Input
% X             -> Grid of Tx Symbols
% Y             -> Grid of Rx Symbols


nREs = size(Y,1);
nSymbols = size(Y, 2);
nLayers = size(Y, 3);

nRBs = round(nREs / 12);


H = zeros(nRBs*4, 10);
H_grid = zeros(nREs, 10);

% Error checking
if (nSymbols ~= 140)
    error('Invalid number of OFDM symbols in the grid');
end

iReCSRsEven = 1:6:nREs;
iReCSRsOdd = 4:6:nREs;

iReCSRs = union(iReCSRsEven, iReCSRsOdd);
%% Estimate channel response
for iSubframe = 1:10 % For each subframe

    % Absolute Symbol index within the frame
    iEvenSymbol = (iSubframe - 1)*14 + [1 8];
    iOddSymbol = (iSubframe - 1)*14 + [5 12];

    H(1:2:end, iSubframe) = sum(Y(iReCSRsEven, iEvenSymbol) ...
        ./ X(iReCSRsEven, iEvenSymbol), 2) / 2;
    H(2:2:end, iSubframe) = sum(Y(iReCSRsOdd, iOddSymbol) ...
        ./ X(iReCSRsOdd, iOddSymbol), 2) / 2;

end

%% Interpolate

H_grid(iReCSRs, :) = H;

L_mov = 4;
h_mov = (1/L_mov)*ones(L_mov, 1);

%% Obtain FEQ

FEQ = 1./H_grid;
end
