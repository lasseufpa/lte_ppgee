function csr = CSRgenerator(iSubframe, numTx, nRBs)
% LTE Cell-Specific Reference signal generation. Section 6.10.1 of 3GPP TS
% 36.211 v10.0.0. Generate the whole set per OFDM symbol, for 2 OFDM
% symbols per slot, for 2 slots per subframe, per antenna port (numTx).
% This fcn accounts for the per antenna port sequence generation, while the
% actual mapping to resource elements is done in the Resource mapper.
%   Input:
% iSubframe --> Subframe index
% numTx     --> Number of Tx antennas
% nRBs      --> Number of Resource Blocks
%
%   Output:
% csr       --> CSR Symbols with dim (subcarrier, symbol, slot, antenna)
% Notes:
%   1)  The "subcarrier" dimension in csr refers to the subcarriers that
%       carry CSR only. It is twice the number of Resource Blocks because
%       for each OFDM symbol, each RB has two CSR symbols.
%
%   2) The "symbol" and "slot" numer of dimensions are both always 2,
%       because there are 2 OFDM symbols with CSR per slot and 2 slots per
%       subframe. The product of the two dimensions gives the 4 OFDM
%       symbols carrying CSR.
%
% %#codegen
persistent hSeqGen;
persistent hInt2Bit;
% Assumed parameters
NcellID = 0; % One of possible 504 values
Ncp = 1; % for normal CP, or 0 for Extended CP
csr = complex(zeros(nRBs*2, 2, 2, numTx));
l = [0; 4]; % OFDM symbol idx in a slot for common first antenna port
% Buffer for sequence per OFDM symbol
seq = zeros(size(csr,1)*2, 1); % *2 for complex outputs
if isempty(hSeqGen)
    hSeqGen = comm.GoldSequence(...
        'FirstPolynomial',[1 zeros(1, 27) 1 0 0 1],...
        'FirstInitialConditions', [zeros(1, 30) 1], ...
        'SecondPolynomial', [1 zeros(1, 27) 1 1 1 1],...
        'SecondInitialConditionsSource', 'Input port',...
        'Shift', 1600,...
        'SamplesPerFrame', length(seq));
    hInt2Bit = comm.IntegerToBit('BitsPerInteger', 31);
end
% Generate the common first antenna port sequences
for i = 1:2 % slot wise
    for lIdx = 1:2 % symbol wise
        c_init = (2^10)* (7*((iSubframe+i-1)+1)+l(lIdx)+1)*(2*NcellID+1)...
            + 2*NcellID + Ncp;
        % Convert to binary vector
        iniStates = step(hInt2Bit, c_init);
        % Scrambling sequence - as per Section 7.2, 36.211
        seq = step(hSeqGen, iniStates);
        % Store the common first antenna port sequences
        csr(:, lIdx, i, 1) = ...
            (1/sqrt(2))*complex(1-2.*seq(1:2:end), 1-2.*seq(2:2:end));
    end
end
% Copy the duplicate set for second antenna port, if exists
if (numTx>1)
    csr(:, :, :, 2) = csr(:, :, :, 1);
end
% Also generate the sequence for l=1 index for numTx = 4
if (numTx>2)
    for i = 1:2 % slot wise %l=1
        c_init = (2^10)*(7*((iSubframe+i-1)+1)+1+1)*(2*NcellID+1) + ...
            2*NcellID + Ncp;
        % Convert to binary vector
        iniStates = step(hInt2Bit, c_init);
        % Scrambling sequence - as per Section 7.2, 36.211
        seq = step(hSeqGen, iniStates);
        % Store the third antenna port sequences
        csr(:, 1, i, 3) = ...
            (1/sqrt(2))*complex(1-2.*seq(1:2:end), 1-2.*seq(2:2:end));
    end
    % Copy the duplicate set for fourth antenna port
    csr(:, 1, :, 4) = csr(:, 1, :, 3);
end