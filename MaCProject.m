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
Nbps= 2;                    % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6)
CutoffFreq= 1000000;        % CutOff Frequency of the Nyquist Filter
RollOff= 0.3;               % Roll-Off Factor
OSF= 4;                     % Oversampling Factor
Tsymb= 1/(2*CutoffFreq);    % Symbol Period
BitRate= Nbps/Tsymb;        % Bit Rate
SymRate= 1/Tsymb;           % Symbol Rate
Fs = OSF*SymRate;           % Sampling Frequency
N = 101;                    % Number of taps (ODD ONLY)
EbN0 = 100;                   % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)

%=============================================%
% Bit Generation
%------------------------

bits_tx = (randi(2,1,Nb)-1)';               % bits_tx = Binary sequence

% Mapping
%------------------------
if Nbps>1
    In_tx = mapping(bits_tx,Nbps,'qam');         % In = Symbols sequence at transmitter
else
    In_tx = mapping(bits_tx,Nbps,'pam');         
end

figure;
plot(In_tx,'r*');
hold off;

% Upsampling
%-----------------------
signal = zeros(Nb/Nbps*OSF,1);

for i = 1:Nb/Nbps
    signal(1+OSF*(i-1))=In_tx(i);
    for j = 2:OSF
        signal(j+OSF*(i-1))=0;
    end
end

% RRC Nyquist Filter TX
%-------------------------

[h_RRC,H_RRC] =  RRC(Fs,Tsymb,N,RollOff);
signal_tx = conv(signal,h_RRC);

figure("Name","Symbol sequence at transmitter");
plot(signal_tx,"r*")
hold off;

% Noise
%------------------

SignalEnergy = (trapz(abs(In_tx).^2))*(1/Fs);
Eb = SignalEnergy/(2*Nb);

N0 = Eb/(10.^(EbN0/10));
NoisePower = 2*N0*Fs;

noise = (sqrt(NoisePower/2)*(randn(1,length(signal_tx))+1i*randn(1,length(signal_tx))))';

signal_rx = signal_tx + noise;

figure("Name","Noised RX signal");
plot(signal_rx,"r*")
hold off;

% RRC Nyquist Filter RX
%-------------------------

filtered_signal_rx = conv(signal_rx,h_RRC);
filtered_signal_rx = filtered_signal_rx(N:end-(N-1));

figure("Name","Filtered RX signal");
plot(filtered_signal_rx,"r*")
hold off;

% Downsampling
%-------------

downsampled_signal = zeros(Nb/Nbps,1);

for i = 1:Nb/Nbps
    downsampled_signal(i)=sum(filtered_signal_rx(1+OSF*(i-1):i*OSF))/OSF;
end


% Demapping
%-----------

if Nbps>1
    bits_rx = demapping(downsampled_signal,Nbps,"qam");
else
    bits_rx = demapping(downsampled_signal,Nbps,"pam");
end

% BER
%-----

ber = 0; 
for k=1:length(y)
    if(bits_rx(k) ~= b(k))
        ber = ber+1;
    end
end
ber = ber/(Nb/Nbps)
% % %-------------------------------------%
