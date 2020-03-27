%-------------------------------------%
%    Modulation and Coding Project    %
%-------------------------------------%
%   Authors : Theo LEPOUTTE           %
%             John ANTOUN             %
%                                     %
%   Date : March 16, 2020             %
%-------------------------------------%
clc; clear all; close all;
%------Parameters------%

Nb= 1200;                   % Number of bits    
Nbps= 4;                    % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6)
CutoffFreq= 1000000;        % CutOff Frequency of the Nyquist Filter
RollOff= 0.3;               % Roll-Off Factor
OSF= 4;                     % Oversampling Factor
Tsymb= 1/(2*CutoffFreq);    % Symbol Period
BitRate= Nbps/Tsymb;        % Bit Rate
SymRate= 1/Tsymb;           % Symbol Rate
Fs = OSF*SymRate;           % Sampling Frequency
N = 101;                    % Number of taps (ODD ONLY)
EbN0 = 10;                   % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)

%=============================================%
% Bit Generation
%------------------------

bits_tx = randi(2,1,Nb)-1;               % bits_tx = Binary sequence

% Mapping
%------------------------
if Nbps>1
    signal_tx = mapping(bits_tx.',Nbps,'qam')';         % In = Symbols sequence at transmitter
else
    signal_tx = mapping(bits_tx.',Nbps,'pam')';         
end

% Upsampling
%-----------------------

upsampled_signal = zeros(1,Nb/Nbps*OSF);
for i = 1:Nb/Nbps
    upsampled_signal(1+OSF*(i-1))=signal_tx(i);
    for j = 2:OSF
        upsampled_signal(j+OSF*(i-1))=0;
    end
end

% RRC Nyquist Filter TX
%-------------------------

[h_RRC,H_RRC] =  RRC(Fs,Tsymb,N,RollOff);
filtered_signal_tx = conv(upsampled_signal,h_RRC);

figure("Name","TX signal")
subplot(1,2,1)
plot(upsampled_signal,'r*')
title("Upsampled TX signal")
subplot(1,2,2)
plot(filtered_signal_tx,"r*")
title("Filtered TX signal")

% Noise
%------------------

SignalEnergy = (trapz(abs(filtered_signal_tx).^2))*(1/Fs);
Eb = SignalEnergy/(2*Nb);

N0 = Eb/(10^(EbN0/10));
NoisePower = 2*N0*Fs;

noise = sqrt(NoisePower/2)*(randn(1,length(filtered_signal_tx))+1i*randn(1,length(filtered_signal_tx)));

signal_rx = filtered_signal_tx + noise;

figure("Name","RX signal");
subplot(1,2,1);
plot(signal_rx,"r*")
title("Noised RX signal");

% RRC Nyquist Filter RX
%-------------------------

filtered_signal_rx = conv(signal_rx,fliplr(h_RRC));
cropped_filtered_signal_rx = filtered_signal_rx(N:end-(N-1));

subplot(1,2,2);
plot(cropped_filtered_signal_rx,"r*")
title("Filtered RX signal");

% Downsampling
%-------------

downsampled_signal = zeros(1,Nb/Nbps);
for i = 1:Nb/Nbps
    downsampled_signal(i)=cropped_filtered_signal_rx(1+OSF*(i-1));
end

% Demapping
%-----------

bits_rx = zeros(1,Nb);
if Nbps>1
    bits_rx = demapping(downsampled_signal',Nbps,"qam");
else
    bits_rx = demapping(real(downsampled_signal'),Nbps,"pam");
end

% BER
%-----

BER = 0; 
for i=1:Nb
    if(bits_rx(i) ~= bits_tx(i))
        BER = BER+1;
    end
end
BER = BER/Nb
 
