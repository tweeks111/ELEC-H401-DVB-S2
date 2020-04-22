%-------------------------------------%
%    Modulation and Coding Project    %
%-------------------------------------%
%   Authors : Theo LEPOUTTE           %
%             John ANTOUN             %
%                                     %
%   Date : March 16, 2020             %
%-------------------------------------%
clc;clear;close all;
addpath('../Part 1 - Communication Chain');
%------Parameters------%
Nb= 600;                                        % Number of bits  
Nbps= 4;                                        % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6) -> vector to compare 
CutoffFreq= 1e6;                                % CutOff Frequency of the Nyquist Filter
RollOff= 0.3;                                   % Roll-Off Factor
M= 4;                                           % Upsampling Factor
N = 23;                                         % Number of taps (ODD ONLY)
EbN0 = 10;                                      % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
Tsymb= 1/(2*CutoffFreq);                        % Symbol Period
SymRate= 1/Tsymb;                               % Symbol Rate
Fs = SymRate*M;                                 % Sampling Frequency

%=============================================%
% Bit Generation
%------------------------

bits_tx = randi(2,1,Nb)-1;               % bits_tx = Binary sequence

% Mapping
%------------------------
if Nbps(nbps)>1
    signal_tx = mapping(bits_tx.',Nbps(nbps),'qam').';         % Symbols sequence at transmitter
else
    signal_tx = mapping(bits_tx.',Nbps(nbps),'pam').';         
end

% Upsampling
%-----------------------

