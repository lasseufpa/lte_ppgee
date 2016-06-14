function [ bn, bn_bar ] = bitloadResourceBlock( b )
% Generates an OFDM loading vector based on the bit loads adopted for each
% resource block.
%
% Input
% b     -> Bit loading of resource blocks
%
% Output
%
% bn        -> Bits per complex subchannel
% bn_bar    -> Bits per real dimension of each subchannel

%% Constants

nSubchannelsPerRb = 12;
nRealDimensionsPerOfdmSubchannel = 2;

%% Processing

% Vector of bits per OFDM subchannel (bits per complex dimension)
bn = repmat(b, nSubchannelsPerRb, 1);
bn = bn(:);

% Vector of bits per real dimension of the OFDM subchannels
%
% Note each subchannel is two-dimensional and all subchannels in a RB are
% loaded with the same number of bits.
bn_bar = bn / nRealDimensionsPerOfdmSubchannel;

end

