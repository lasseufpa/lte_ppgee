function [ params ] = lteOfdmParameters( bw_choice )
% Returns OFDM parameters for each LTE BW choice
% Note: the FFT size is not standardized, but often the bandwidths below
% are associated with the given FFT sizes.

switch (bw_choice)
    case 1
        % LTE 1.4 MHz
        nRBs = 6;
        W    = 1.4e6;
        fs   = 1.92e6;
        N    = 128;
    case 2
        % LTE 3 MHz
        nRBs = 15;
        W    = 3e6;
        fs   = 3.84e6;
        N    = 256;
    case 3
        % LTE 5 MHz
        nRBs = 25;
        W    = 5e6;
        fs   = 7.68e6;
        N    = 512;
    case 4
        % LTE 10 MHz
        nRBs = 50;
        W    = 10e6;
        fs   = 15.36e6;
        N    = 1024;
    case 5
        % LTE 15 MHz
        nRBs = 75;
        W    = 15e6;
        fs   = 23.04e6;
        N    = 1536;
    case 6
        % LTE 20 MHz
        nRBs = 100;
        W    = 20e6;
        fs   = 30.72e6;
        N    = 2048;
end

%% Cp Length
normalCp = 4.68e-6;
biggerCp = 5.2e-6;

nu     = round(normalCp * fs); % Normal CP
nu_ext = round(biggerCp * fs); % Extended CP

%% Output parameters
params.nUsedSubcarriers = nRBs * 12; % Used subcarriers
params.BW               = W;         % Bandwidth
params.fs               = fs;        % Sampling Frequency
params.N                = N;         % FFT size
params.nu               = nu;        % Normal CP
params.nu_ext           = nu_ext;    % Extended CP

end

