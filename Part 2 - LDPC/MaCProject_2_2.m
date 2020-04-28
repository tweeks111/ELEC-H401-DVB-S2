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

%%
% Bit Generation
%------------------------

bits_tx = randi(2,1,Nb)-1;               % bits_tx = Binary sequence

%%
% LDPC
%----------------

H0 = makeLDPC(128, 256,0,1,3);
[codedbits, H] = makeParityChk(bits_tx, H0, 0);


%%
% Mapping
%------------------------
if Nbps>1
    signal_tx = mapping(codedbits.',Nbps,'qam').';         % Symbols sequence at transmitter
else
    signal_tx = mapping(codedbits.',Nbps,'pam').';         
end


%%
% Upsampling
%-----------------------
upsampled_signal = zeros(1,Nb/Nbps*M);
for i = 1:Nb/Nbps
    upsampled_signal(1+M*(i-1))=signal_tx(i);
    for j = 2:M
        upsampled_signal(j+M*(i-1))=0;
    end
end

%%
% RRC Nyquist Filter TX
%-------------------------
[h_RRC,H_RRC] =  RRC(Fs,Tsymb,N,RollOff,Nbps,AverageNb,M);
filtered_signal_tx = conv(upsampled_signal,h_RRC);


%%
% Noise
%-----------------

SignalEnergy = (trapz(abs(filtered_signal_tx).^2))*(1/Fs);
Eb = SignalEnergy/(2*Nb);

N0 = Eb./(10.^(EbN0/10));
NoisePower = 2*N0*Fs;

noise = zeros(length(EbN0),Nb/Nbps*M+N-1);
signal_rx = zeros(length(EbN0),Nb/Nbps*M+N-1);
for j = 1:length(EbN0)
    noise(j,:) = sqrt(NoisePower(j)/2).*(randn(1,Nb/Nbps*M+N-1)+1i*randn(1,Nb/Nbps*M+N-1));
    signal_rx(j,:) = filtered_signal_tx + noise(j,:);
end


%%
% RRC Nyquist Filter RX
%-------------------------

filtered_signal_rx = zeros(length(EbN0),Nb/Nbps*M+2*(N-1));
cropped_filtered_signal_rx = zeros(length(EbN0),Nb/Nbps*M);
for i =1:length(EbN0)
    filtered_signal_rx(i,:) = conv(signal_rx(i,:),fliplr(h_RRC));
    cropped_filtered_signal_rx(i,:) = filtered_signal_rx(i,N:end-(N-1));
end

%%
% Downsampling
%-------------

downsampled_signal = zeros(length(EbN0),Nb/Nbps);
for j = 1:length(EbN0)
    for i = 1:Nb/Nbps
        downsampled_signal(j,i)=cropped_filtered_signal_rx(j,1+M*(i-1));
    end
end

%%
% Demapping
%-----------

codedbits_rx = zeros(length(EbN0),Nb);
for j = 1:length(EbN0)
    if Nbps>1
        codedbits_rx(j,:) = demapping(downsampled_signal(j,:).',Nbps,"qam");
    else
        codedbits_rx(j,:) = demapping(real(downsampled_signal(j,:).'),Nbps,"pam");
    end
end

%%
% Soft Decoding
%----------------

variance=NoisePower/2;
bits_rx=softDecoding(codedbits_rx,H,variance);




