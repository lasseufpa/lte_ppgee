function [ FEQ ] = feqFromCsrSymbols(X, Y)
% Estimate FEQ from CSR Symbols
%
% This function extracts the CSRs from the Rx grid and compares to the CSRs
% that should be obtained (as in the Tx grid) to the obtain the FEQ and
% sparse points in the grid. It, then, interpolates the FEQ to obtain the
% equalizer gain in all Resource Elements. In the end, one FEQ matrix is
% output with number of rows equivalent to the number of resource elements
% and number of columns equivalent to the number of symbols in a frame,
% even though the same FEQ is applied to all symbols in a subframe. The
% third dimension of the output FEQ matrix corresponds to the layer.
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

% Type of Interpolation
interpType = 1; % 0 - MovingAverage; 1 - Linear

% Constants
nSubframes = 10;

% Infer Parameters
nREs     = size(Y, 1);
nSymbols = size(Y, 2);
nLayers  = size(Y, 3);
nRBs    = round(nREs / 12);

% Preallocate
H_est   = zeros(nRBs*4, nSubframes, nLayers); % Only the CSR subcarriers
H_grid  = zeros(nREs, nSubframes, nLayers);   % All subcarriers
FEQ = zeros(nREs, 140, nLayers);          % FEQ Output

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

for iLayer = 1:nLayers
    for iSubframe = 1:10 % For each subframe

        % Absolute Symbol index within the frame
        iEvenSymbol = (iSubframe - 1)*14 + [1 8];
        iOddSymbol  = (iSubframe - 1)*14 + [5 12];

        % Channel tap estimation:
        H_est(1:2:end, iSubframe, iLayer) = ...
            sum(Y(iReCSRsEven, iEvenSymbol, iLayer) ...
            ./ X(iReCSRsEven, iEvenSymbol, iLayer), 2) / 2;
        H_est(2:2:end, iSubframe, iLayer) = ...
            sum(Y(iReCSRsOdd, iOddSymbol, iLayer) ...
            ./ X(iReCSRsOdd, iOddSymbol, iLayer), 2) / 2;

    end
end

%% Interpolate
% There are only 4 Resource Elements per Resource Block with CSRs. Thus,
% the number of channel taps that were estimated for each subframe and RB
% was 4. We need now to fill all the 12 REs of each RB with estimated taps.
% This is done by interpolation.

% Fill the matrix that contains the channel estimation at all subcarriers,
% with the CSR estimations:
H_grid(iReCSRs, :) = H_est;

% Effectively interpolate:
switch (interpType)
    case 0
        % Replicate the taps for each three neighbor subcarriers
        L_mov = 3;
        h_mov = ones(1, L_mov);
        H_grid = filter(h_mov, 1, H_grid);
    case 1
        % Apply linear interpolation
        % This essentialy draws a line between every two CSR estimations
        % and computes the intermediate points in such lines.
        H_grid = interp1q([iReCSRs.'; iReCSRs(end)+3],...
            [H_grid(iReCSRs, :); H_grid(iReCSRs(end), :)],...
            (1:nREs).');
end

%% Obtain FEQ
% Output a matrix of 140 FEQs

% FEQ for each subframe:
FEQ_per_subframe = 1./H_grid;

% Replicate the FEQ for each subframe such that a N x 140 matrix is
% generated:
for iSubframe = 1:nSubframes
    iSymbol = (iSubframe - 1)*14 + (1:14);
    FEQ(:, iSymbol) = repmat(FEQ_per_subframe(:,iSubframe), [1 14]);
end

end
