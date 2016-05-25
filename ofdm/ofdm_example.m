%% Simple OFDM Example
clearvars, clc

%% Debug levels
debug               = 1;  % Enable debug information
debug_constellation = 0;  % Debug a certain subchannel constellation
debug_tone          = 16; % Tone whose constellation is debugged
debug_Pe            = 1;  % Debug error probabilities
debug_tx_energy     = 0;  % Debug transmit energy
debug_snr           = 0;  % Debug SNRs

%% Parameters
alpha      = 1;         % Increase FFT size by this factor preserving Fs
% Note: this is useful to evaluate the DMT performance as N -> infty
Px         = 1e-3;      % Transmit Power (W)
N0_over_2  = 1e-14;     % Noise PSD (W/Hz/dim) and variance per dimension
N          = 1024;      % FFT size and the number of used real dimensions
nu         = 80;        % Cyclic Prefix Length
delta_f    = 15e3;      % Subchannel bandwidth
nSymbols   = 1e3;       % Number of DMT symbols per transmission iteration
subchanOrder = 16;       % Modulation order adopted for all subchannels
%   Note: LTE allows 4-QAM, 16-QAM and 64-QAM
% Monte-Carlo Parameters
maxNumErrs   = 1e5;
maxNumDmtSym = 1e12;

%% Derived computations:
N         = N*alpha;
delta_f   = delta_f/alpha;
Fs        = N * delta_f;
nDim      = 2*(N + nu);     % Total number of real dimensions per OFDM
% symbol
Ts        = 1 / Fs;
Tsym      = (N + nu) * Ts;  % Symbol Period
Rsym      = 1 / Tsym;       % DMT Symbol rate (real dimensions per sec)
Ex        = Px * Tsym;      % Average DMT symbol energy
Ex_bar    = Ex / nDim;      % Energy per real dimension
% Note the degrees of freedom in the cyclic prefix are taken into account
% when computing the number of dimensions in which the transmit energy is
% distributed.

%% System Objects

EVM = comm.EVM;

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

%% Bit and Energy loading per subchannel

% Vector of bits per dimension in each subchannel (not each subchannel is
% two-dimensional):
bn_bar = log2(subchanOrder) * ones(N,1);
% Vector of bits per subchannel
bn = 2*bn_bar;

% Energy per dimension in each subchannel (flat energy load):
En_bar = Ex_bar * ones(N,1);
% Energy per two-dimensional subchannel:
En = 2*En_bar;

%% SNR per subchannel

% Note: the SNR in each subchannel is given by:
SNR_n = En_bar .* gn;

if (debug_snr)
    figure
    plot(10*log10(SNR_n))
    hold on
    plot(10*log10(SNRmfb) * ones(N, 1), '--r')
    xlabel('Tone (n)')
    ylabel('SNR_n (db)')
    set(gca,'XLim',[1 N]);
    drawnow 
end

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

% Number of bits per two-dimension (complex sample)
rho = (1/(nDim/2))*(sum(bn));
fprintf('Spectral efficiency:    \t %g bits/2-dimensions\n', rho)
% For gap=0 and N->+infty, this should be the channel capacity per real
% dimension.

% Corresponding multi-channel SNR:
SNRdmt    = gap*(2^rho - 1);
SNRdmt_db = 10*log10(SNRdmt);
% SNRdmt is the SNR that reflects the achieved bit-rate considering the
% given gap to capacity. For a gap of 0, SNRdmt approaches the channel
% capacity as N (FFT size) goes to infinity.
fprintf('Multi-channel SNR (SNRdmt):\t %g dB\n', SNRdmt_db)

% Bit-rate
Rbit = rho * Fs; % (bits/2-Dim) * (2-Dim/sec) = bits/sec
% This should be equivalent to "Fs * log2(1 + SNRdmt/gap)"
% Capacity
c = Fs * log2(1 + SNRdmt);
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

fprintf('Expected EVM:\t                %g %%rms\n', 100 * expectedEVM);

%% Modulators

% Modulation orders per subchannel
modOrder = 2.^bn;
% Unique modulation orders of the two-dimensional subchannels:
twoDim_const_orders = unique(modOrder(modOrder~=1));

%Preallocate modems
modulator = cell(length(modOrder), 1);
demodulator = cell(length(modOrder), 1);

% Configure 2-dimensional modems for each distinct bit loading:
for i = 1:length(twoDim_const_orders)
    M = twoDim_const_orders(i);
    
    modulator{i} = modem.qammod('M', M, 'SymbolOrder', 'Gray');
    demodulator{i} = modem.qamdemod('M', M, 'SymbolOrder', 'Gray');
end

%% Look-up table for each subchannel indicating the corresponding modem

modem_n = zeros(N, 1);

for k = 1:N
    iModem = find(twoDim_const_orders == modOrder(k));
    if (iModem)
        modem_n(k) = iModem;
    end
end

%% Energy loading (constellation scaling factors)

Scale_n = zeros(N, 1);

for k = 1:N
    Scale_n(k) = modnorm(...
        modulator{modem_n(k)}.constellation,...
        'avpow', En(k));
end

%% Monte-carlo

fprintf('\n---------------------- Monte Carlo --------------------- \n\n');

% Preallocate
X          = zeros(N, nSymbols);
tx_symbols = zeros(N, nSymbols);
rx_symbols = zeros(N, nSymbols);
sym_err_n  = zeros(N, 1);

numErrs = 0; numDmtSym = 0;

% Sys Objects
BitError = comm.ErrorRate;

%% Printing Header

if (debug && debug_tx_energy)
    fprintf('Tx Energy   |\t');
    fprintf('Nominal     |\t');
end
fprintf('SNR (Time)  |\t');
fprintf('RMS EVM     |\t');
fprintf('Pe_bar      |\t');
fprintf('nErrors     |\t');
fprintf('nDMTSymbols |\t\n');

%% Iterative Transmissions

iTransmission = 0;

while ((numErrs < maxNumErrs) && (numDmtSym < maxNumDmtSym))
    iTransmission = iTransmission + 1;
    
    % Random Symbol generation
    for k = 1:N
        tx_symbols(k, :) = randi(modOrder(k), 1, nSymbols) - 1;
    end
    
    %% Constellation Encoding
    for k = 1:N
        if (modem_n(k) > 0)
            X(k, :) = Scale_n(k) * ...
                modulator{modem_n(k)}.modulate(tx_symbols(k, :));
        end
    end
    
    %% Modulation
    
    x = sqrt(N) * ifft(X, N); % Orthonormal IFFT
    
    %% Cyclic extension
    
    x_ext = [x(N-nu+1:N, :); x];
    
    %% Parallel to serial
    
    u = x_ext(:);
    
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
            tx_total_energy / nSymbols);
        fprintf('%12g|\t',Ex);
    end
    
    %% Channel    
    y = conv(u, p);    
    
    %% Noise
    
    nn = sqrt(N0_over_2) * (randn(length(y),1) + 1j*randn(length(y),1));
    
    % Add noise
    y = y + nn;
    
    % SNR in time-domain:
    SNR_time = 10*log10(mean(abs(y).^2) / mean(abs(nn).^2));
    
    fprintf('%12g|\t', SNR_time);
    
    %% Synchronization
    % Note: synchronization introduces a phase shift that should be taken
    % into account in the FEQ.
    
    nRxSamples = (N+nu)*nSymbols;
    y_sync     = y((n0 + 1):(n0 + nRxSamples));
    
    %% Serial to Parallel
    
    y_sliced = reshape(y_sync, N + nu, nSymbols);
    
    %% Extension removal
    
    y_no_ext = y_sliced(nu + 1:end, :);
    
    %% Regular Demodulation (without decision feedback)
    
    % FFT
    Y = (1/sqrt(N)) * fft(y_no_ext, N); % Orthonormal FFT
    
    % FEQ - One-tap Frequency Equalizer
    Z = diag(FEQ) * Y;
    
    %% EVM
    RMSEVM = step(EVM, X(:), Z(:));
    
    fprintf('%12g|\t', RMSEVM);
    
    %% Constellation decoding (decision)
    
    for k = 1:N
        if (modem_n(k) > 0)
            rx_symbols(k, :) = demodulator{modem_n(k)}.demodulate(...
                (1/Scale_n(k)) * Z(k, :));
        end
    end
    
    % Symbol error count
    sym_err_n = sym_err_n + symerr(tx_symbols, rx_symbols, 'row-wise');
    % Symbol error rate per subchannel
    ser_n     = sym_err_n / (iTransmission * nSymbols);
    % Per-dimensional symbol error rate per subchannel
    ser_n_bar = ser_n / 2;
    
    % Preliminary results
    numErrs   = sum(sym_err_n);
    numDmtSym = iTransmission * nSymbols;
    
    fprintf('%12g|\t', mean(ser_n_bar));
    fprintf('%12g|\t', numErrs);
    fprintf('%12g|\n', numDmtSym);
    
    
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
