%% LTE OFDM Simulator
clearvars, clc

%% Debug levels
debug               = 1;  % Enable debug information
debug_constellation = 0;  % Debug a certain subchannel constellation
debug_tone          = 16; % Tone whose constellation is debugged
debug_Pe            = 1;  % Debug error probabilities
debug_tx_energy     = 0;  % Debug transmit energy
debug_snr           = 0;  % Debug SNRs
debug_evm           = 0;  % Debug EVMs
debug_psd           = 0;  % Debug the Rx signal PSD

%% Parameters
lteChoice  = 5;
% 1 -> LTE 1.4
% 2 -> LTE 3
% 3 -> LTE 5
% 4 -> LTE 10
% 5 -> LTE 15
% 6 -> LTE 20
Px         = 1e-3;      % Transmit Power (W)
N0_over_2  = 1e-14;     % Noise PSD (W/Hz/dim) and variance per dimension
delta_f    = 15e3;      % Subchannel bandwidth
nSymbols   = 1e3;       % Number of OFDM symbols per transmission iteration
subchanOrder = 16;      % Modulation order adopted for all subchannels
nLayers      = 1;       % Number of transmit layers
%   Note: LTE allows 4-QAM, 16-QAM and 64-QAM
% Monte-Carlo Parameters
maxNumErrs   = 1e5;
maxNumOfdmSym = 1e12;

%% Derived computations:

% Get OFDM parameters that depend on the LTE choice
[ lte ] = lteOfdmParameters( lteChoice );
N       = lte.N;
N_used  = lte.nUsedSubcarriers;
nRBs    = lte.nRBs;
nu      = lte.nu;
Fs      = lte.fs;

% Number of OFDM Symbols per LTE Frame
nSymbolsPerFrame   = 140;

% Index of the subchannels that should be loaded:
used_tones = [ 2:(N_used/2), (N - N_used/2):N].';
% Note DC is not used and only the bands adjcent to DC are used. Recall the
% FFT order is such that the first N/2 + 1 tones represent the "positive
% spectrum" (0 to pi) and the last N/2 -1 tones represent the negative half
% of the baseband-centered spectrum (pi to 2*pi).


nDim      = 2*(N + nu);     % Total number of real dimensions per OFDM
% symbol
Ts        = 1 / Fs;
Tsym      = (N + nu) * Ts;  % Symbol Period
Rsym      = 1 / Tsym;       % OFDM Symbol rate (real dimensions per sec)
Ex        = Px * Tsym;      % Average OFDM symbol energy
Ex_bar    = Ex / nDim;      % Energy per real dimension
% Note the degrees of freedom in the cyclic prefix are taken into account
% when computing the number of dimensions in which the transmit energy is
% distributed.

%% System Objects

EVM = comm.EVM;
EVM.AveragingDimensions = [1 1];
SpecAnalyzer = dsp.SpectrumAnalyzer('SampleRate',lte.fs);

%% Pulse Response

% "Ideal" flat-fading Channel
p = [1];

% Pulse frequency response
H = fft(p, N);
H = H(:);

% Pulse response length
Lh = length(p);

% Matched-filter Bound
SNRmfb = (Ex_bar * norm(p).^2) / N0_over_2;
fprintf('SNRmfb:    \t %g dB\n\n', 10*log10(SNRmfb))

%% 1-tap Frequency Equalizer

% Frequency domain response
H_freq = fft(p, N);

% Cursor
[~, iMax] = find(abs(p)>1e-5, 1, 'first');
n0 = iMax - 1;

% Corresponding phase shift due to cursor
phaseShift = exp(1j*2*pi*(n0/N)*(0:N-1));

% FEQ
FEQ = 1 ./ (H_freq .* phaseShift);

%% SNR for unitary energy transmit symbol
% This is also called the "channel gain to noise ratio"

gn = (abs(H).^2) / N0_over_2;

%% Bit load per Resource Block and per Subchanne;

% Vector of bits per subchannel in each RB (bits per complex dimension)
b_nRb = log2(subchanOrder) * ones(nRBs, 1);

% Obtain the bit loading vectors:
% bn        -> Bits per complex subchannel
% bn_bar    -> Bits per real dimension of each subchannel
[bn, bn_bar] = bitloadResourceBlock( b_nRb );

%% Energy load per subchannel

% 1) What is the budget of energy available for each OFDM symbol? Ex
% 2) How many subchannels are effectively loaded with energy? N_used
% Thus, assuming all subchannels are loaded with an equal amount of energy,
% the energy load per subchannel is:
En_per_used_subchannel = (Ex / N_used);
% Then, since each subchannel is two-dimensional in OFDM, the energy per
% dimension in each n-th subchannel becomes:
En_bar = (En_per_used_subchannel/2) * ones(N_used,1);
% Finally, since the CP uses part of the transmit energy, the energy per
% subchannel has to discount the CP repetition so that the transmit energy
% budget is obeyed:
En_bar = (N/(N+nu)) * En_bar;
% Energy per two-dimensional subchannel:
En = 2*En_bar;

%% SNR per subchannel

% Note: the SNR in each subchannel is given by:
SNR_n = En_bar .* gn(used_tones);

%% Gap to capacity, Multi-channel SNR and Channel Capacity

% With this per-dimensional SNR, we can compute the gap to capacity in each
% subchannel. Recall the gap to capacity, given by:
%
%   gap = SNRnorm_n = SNR_n/(2^rho - 1),
%
% where rho is the spectral efficiency in bits per complex dimension, which
% is equivalent to 2*bn_bar.
gap_n     = SNR_n ./(2.^(2*bn_bar) - 1);
gap_db_n  = 10*log10(gap_n); % in dB

gap = mean(gap_n);
gap_db = 10*log10(gap);

fprintf('Gap to capacity:        \t %g db\n', gap_db)

% Total number of bits per two-dimension (complex sample), i.e. the
% spectral effiency:
rho = (1/(nDim/2))*(sum(bn));
fprintf('Spectral efficiency:    \t %g bits/2-dimensions\n', rho)
% For gap -> 0 and N -> +infty, this should be the channel capacity per
% real dimension.

% Corresponding multi-channel SNR:
SNRofdm    = gap*(2^rho - 1);
SNRofdm_db = 10*log10(SNRofdm);
% SNRofdm is the SNR that reflects the achieved bit-rate considering the
% given gap to capacity. For a gap of 0, SNRofdm approaches the channel
% capacity as N (FFT size) goes to infinity.
fprintf('Multi-channel SNR (SNRofdm):\t %g dB\n', SNRofdm_db)

% Bit-rate
Rbit = rho * Fs; % (bits/2-Dim) * (2-Dim/sec) = bits/sec
% This should be equivalent to "Fs * log2(1 + SNRofdm/gap)"
% Capacity
c = Fs * log2(1 + SNRofdm);
fprintf('Bit rate:               \t %g mbps\n', Rbit/1e6);
fprintf('Capacity:               \t %g mbps\n', c/1e6);

%% Analysis of the Error Probability per dimension
% Comparison between the water-filling and the discrete loading

fprintf('\n----------------- Error Probabilities ------------------ \n\n');

% QAM-SQ
Pe_bar_n = 2*(1 - 1./(2.^bn_bar)) .* qfunc(sqrt( 3*gap_n) );

fprintf('Approximate NNUB Pe per dimension: %g\n', ...
    mean(Pe_bar_n,'omitnan'));

%% Expected EVM

expectedEVM = 1/sqrt(mean(SNR_n));
% Note, it can be shown that
%
%   mean(SNR_n) = (En_bar * norm(p)^2) / N0_over_2,
%
% where the energy per dimension in each subchannel (En_bar) is constant
% due to the flat energy load. Furthermore, En_bar is given by the
% computation in the "Energy load" section.

fprintf('Expected EVM:\t                %g %%rms\n', 100 * expectedEVM);

%% Modulators

% Modulation orders per resource block
modOrder = 2.^b_nRb;
% Unique modulation orders of the two-dimensional subchannels:
twoDim_const_orders = unique(modOrder(modOrder~=1));

% Preallocate unique modems
modulator = cell(length(twoDim_const_orders), 1);
demodulator = cell(length(twoDim_const_orders), 1);

% Configure two-dimensional modems for each distinct bit load:
for iModem = 1:length(twoDim_const_orders)
    M = twoDim_const_orders(iModem); % Modulation order

    % Modulator/Demodulator Objects:
    modulator{iModem}   = modem.qammod('M', M, 'SymbolOrder', 'Gray');
    demodulator{iModem} = modem.qamdemod('M', M, 'SymbolOrder', 'Gray');
end

%% Look-up table for each RB indicating the corresponding modem

modem_n = zeros(nRBs, 1);

for k = 1:nRBs
    iModem = find(twoDim_const_orders == modOrder(k));
    if (iModem)
        modem_n(k) = iModem;
    end
end

%% Energy loading (constellation scaling factors) per RB

Scale_n = zeros(nRBs, 1);

for k = 1:nRBs
    Scale_n(k) = modnorm(...
        modulator{modem_n(k)}.constellation,...
        'avpow', En(k));
end

%% Monte-carlo

fprintf('\n---------------------- Monte Carlo --------------------- \n\n');

% Preallocate
X          = zeros(N, nSymbolsPerFrame, nLayers);
Z          = zeros(N, nSymbolsPerFrame, nLayers);
tx_symbols = zeros(N_used, nSymbolsPerFrame, nLayers);
rx_symbols = zeros(N_used, nSymbolsPerFrame, nLayers);
sym_err_n  = zeros(N_used, 1);

numErrs = 0; numOfdmSym = 0;

% Sys Objects
BitError = comm.ErrorRate;

%% Printing Header

if (debug && debug_tx_energy)
    fprintf('Tx Energy   |\t');
    fprintf('Nominal     |\t');
end
if (debug && debug_snr)
    fprintf('SNR (Time)  |\t');
end
if (debug && debug_evm)
    fprintf('RMS EVM     |\t');
end
if (debug && debug_Pe)
    fprintf('Pe_bar      |\t');
    fprintf('nErrors     |\t');
    fprintf('OFDMSymbols |\t\n');
end

%% Iterative Transmission of LTE Frames

iTransmission = 0;

while ((numErrs < maxNumErrs) && (numOfdmSym < maxNumOfdmSym))
    iTransmission = iTransmission + 1;

    % Random Symbol generation
    % Iterate over Resource Blocks
    for iRB = 1:nRBs
        iRE = 12*(iRB-1) + 1:12*iRB;
        tx_symbols(iRE, :, :) = ...
            randi(modOrder(iRB), 12, nSymbolsPerFrame, nLayers) - 1;
    end

    %% Constellation Encoding

    % Iterate over Resource Blocks
    for iRB = 1:nRBs
        % Find:
        iRE     = 12*(iRB-1) + 1:12*iRB; % Resource Element Index
        k       = used_tones(iRE);       % Actual subchannel indexes
        tx_data = tx_symbols(iRE, :, :); % Corresponding transmit data
        iModem  = modem_n(iRB);          % Modem used for the RB

        % Modulate and scale the transmit data
        if (modem_n(iRB) > 0)
            for iLayer = 1:nLayers
                X(k, :, iLayer) = Scale_n(iRB) * ...
                    modulator{iModem}.modulate(tx_data(:,:,iLayer));
            end
        end
    end

    %% Modulation

    x = sqrt(N) * ifft(X, N); % Orthonormal IFFT

    %% Cyclic extension

    x_ext = [x(N-nu+1:N, :, :); x];

    %% Parallel to serial
    % Serialize the IFFT samples

    u = reshape(x_ext, (N+nu)*nSymbolsPerFrame, nLayers);
    % Note: should be a matrix with dimensions nSamples x nPorts

    if (debug && debug_tx_energy)
        % Note: "u" should become samples leaving the DAC. In that case,
        % they would repreent coefficients of the sampling theorem's sinc
        % interpolation formula, which is an orthogonal (non-normal)
        % expansion. However, note x comes from an orthonormal expansion,
        % which is the normalized IDFT. Hence, the energy in x at this
        % point is still given simply by:
        tx_total_energy = norm(u).^2;

        % A Ts factor should multiply the norm if u was a vector of samples
        % out of the DAC, but in this case there would be a scaling factor
        % introduced by the DAC anti-imaging LPF. Both would cancel each
        % other.
        fprintf('%12g|\t', ...
            tx_total_energy / nSymbolsPerFrame);
        fprintf('%12g|\t',Ex);
    end

    %% Channel
    y = conv2(u, p);

    %% Noise

    nn = sqrt(N0_over_2) * (randn(length(y), nLayers) + ...
        1j*randn(length(y), nLayers));

    % Add noise
    y = y + nn;

    % SNR in time-domain:
    SNR_time = 10*log10(mean(abs(y).^2) / mean(abs(nn).^2));

    if (debug && debug_snr)
        fprintf('%12g|\t', SNR_time);
    end

    if (debug && debug_psd)
        step(SpecAnalyzer, y);
    end
    %% Timing Synchronization
    % Note: synchronization introduces a phase shift that should be taken
    % into account in the FEQ.

    nRxSamples = (N+nu)*nSymbolsPerFrame;
    y_sync     = y((n0 + 1):(n0 + nRxSamples), :);

    %% Serial to Parallel

    y_sliced = reshape(y_sync, N + nu, nSymbolsPerFrame, nLayers);

    %% Extension removal

    y_no_ext = y_sliced(nu + 1:end, :, :);

    %% Regular Demodulation (without decision feedback)

    % FFT
    Y = (1/sqrt(N)) * fft(y_no_ext, N); % Orthonormal FFT

    % FEQ - One-tap Frequency Equalizer
    for iLayer = 1:nLayers
        Z(:,:,iLayer) = diag(FEQ) * Y(:, :, iLayer);
    end

    %% EVM
    RMSEVM = step(EVM, X(used_tones,:), Z(used_tones, :));

    if (debug && debug_evm)
        fprintf('%12g|\t', RMSEVM);
    end

    %% Constellation decoding (decision)

    for iRB = 1:nRBs
        % Find:
        iRE        = 12*(iRB-1) + 1:12*iRB;     % Resource Element Index
        k          = used_tones(iRE);           % Actual subchannel indexes
        Z_unscaled = (1/Scale_n(iRB)) * Z(k, :, :); % Unscaled Rx Symbols
        iModem     = modem_n(iRB);                  % Modem used for the RB

        % Demodulate unscaled symbols
        if (iModem > 0)
            for iLayer = 1:nLayers
                rx_symbols(iRE, :, iLayer) = ...
                  demodulator{iModem}.demodulate(Z_unscaled(:, :, iLayer));
            end
        end
    end

    %% Error performance evaluation

    % Symbol error count
    for iLayer = 1:nLayers
        sym_err_n = sym_err_n + symerr(tx_symbols(:,:,iLayer), ...
            rx_symbols(:,:,iLayer), 'row-wise');
    end
    % Symbol error rate per subchannel
    ser_n     = sym_err_n / (iTransmission * nSymbolsPerFrame * nLayers);
    % Per-dimensional symbol error rate per subchannel
    ser_n_bar = ser_n / 2;

    % Preliminary results
    numErrs   = sum(sym_err_n);
    numOfdmSym = iTransmission * (nSymbolsPerFrame * nLayers);

    if (debug && debug_Pe)
        fprintf('%12g|\t', mean(ser_n_bar));
        fprintf('%12g|\t', numErrs);
        fprintf('%12g|\n', numOfdmSym);
    end


end

%% Constellation plot for debugging
if (debug && debug_constellation && modem_n(debug_tone) > 0)
    k = debug_tone;
    figure
    plot(Z(k, :), 'o')
    hold on
    plot(Scale_n(k) * ...
        modulator{modem_n(k)}.modulate(0:modOrder(k) - 1),...
        'ro', 'MarkerSize', 8, 'linewidth', 2)
    legend('Rx', 'Tx')
    title(sprintf('Tone: %d', debug_tone));
end

%% Results
fprintf('\n----------------------- Results ------------------------ \n\n');
fprintf('Pe_bar:       \t %g\n', mean(ser_n_bar));

if (debug && debug_Pe)
    figure
    stem(ser_n_bar)
    hold on
    stem(Pe_bar_n, 'g')
    hold on
    title('Results: Pe per dimension')
    xlabel('Subchannel (n)')
    ylabel('$\bar{Pe}(n)$')
    legend('Measured','Theoretical')
    set(gca,'XLim',[1 N/2+1]);
end
