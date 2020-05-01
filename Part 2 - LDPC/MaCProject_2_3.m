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
Nbps= 1;                                        % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6) -> vector to compare 
CutoffFreq= 1e6;                                % CutOff Frequency of the Nyquist Filter
RollOff= 0.3;                                   % Roll-Off Factor
M= 4;                                           % Upsampling Factor
N = 23;                                         % Number of taps (ODD ONLY)
EbN0 = -2:1:14;                                 % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
Tsymb= 1/(2*CutoffFreq);                        % Symbol Period
SymRate= 1/Tsymb;                               % Symbol Rate
Fs = SymRate*M;                                 % Sampling Frequency
BlockSize = 128;
BlockNb=6;
CodeRate = 1/2;
Nb= BlockSize*BlockNb;                         % Number of bits
AverageNb=20;
AverageBER=zeros(1,length(EbN0));
AverageBER_HD=zeros(1,length(EbN0)); 
AverageBER_SD=zeros(1,length(EbN0)); 
H0 = makeLdpc(BlockSize, BlockSize/CodeRate,0,1,3);

for avr = 1:AverageNb
disp(avr);
%%
% Bit Generation
%------------------------

bits_tx = randi(2,1,Nb)-1;               % bits_tx = Binary sequence

%%
% LDPC
%----------------


blocks=reshape(bits_tx,BlockSize,BlockNb);          % on divise le vecteur de bits en matrice de block
[checkbits, H] = makeParityChk(blocks, H0, 0);

blocks=blocks.';
checkbits=checkbits.';

codedbits=horzcat(checkbits,blocks);
codedbits_tx=reshape(codedbits.',[],1);

%%
% Mapping
%------------------------

if Nbps>1
        signal_tx = mapping(codedbits_tx,Nbps,'qam').';         % Symbols sequence at transmitter
        signal_tx_uncoded = mapping(bits_tx.',Nbps,'qam').';
else
        signal_tx = mapping(codedbits_tx,Nbps,'pam').';         % Symbols sequence at transmitter   
        signal_tx_uncoded = mapping(bits_tx.',Nbps,'pam').';
end

%%
% Upsampling
%-----------------

upsampled_signal = zeros(1,length(signal_tx)*M);
upsampled_uncoded_signal = zeros(1,length(signal_tx_uncoded)*M);
for i = 1:length(signal_tx)
    upsampled_signal(1+M*(i-1))=signal_tx(i);
    for j = 2:M
        upsampled_signal(j+M*(i-1))=0;
    end
end
for i = 1:length(signal_tx_uncoded)
    upsampled_uncoded_signal(1+M*(i-1))=signal_tx_uncoded(i);
    for j = 2:M
        upsampled_uncoded_signal(j+M*(i-1))=0;
    end
end
%%
% RRC Nyquist Filter TX
%-------------------------
[h_RRC,H_RRC] =  RRC(Fs,Tsymb,N,RollOff,Nbps,AverageNb,M);
filtered_signal_tx = conv(upsampled_signal,h_RRC);
filtered_signal_tx_uncoded = conv(upsampled_uncoded_signal,h_RRC);

%%
% Noise
%-----------------

SignalEnergy = (trapz(abs(filtered_signal_tx).^2))*(1/Fs);
Eb = SignalEnergy/(2*Nb/CodeRate);

SignalEnergy_uncoded = (trapz(abs(filtered_signal_tx_uncoded).^2))*(1/Fs);
Eb_uncoded = SignalEnergy_uncoded/(2*Nb);

N0 = Eb./(10.^(EbN0/10));
NoisePower = 2*N0*Fs;

N0_uncoded = Eb_uncoded./(10.^(EbN0/10));
NoisePower_uncoded = 2*N0_uncoded*Fs;

noise = zeros(length(EbN0),length(signal_tx)*M+N-1);
signal_rx = zeros(length(EbN0),length(signal_tx)*M+N-1);
noise_uncoded = zeros(length(EbN0),length(signal_tx_uncoded)*M+N-1);
signal_rx_uncoded = zeros(length(EbN0),length(signal_tx_uncoded)*M+N-1);
for j = 1:length(EbN0)
    noise(j,:) = sqrt(NoisePower(j)/2).*(randn(1,length(signal_tx)*M+N-1)+1i*randn(1,length(signal_tx)*M+N-1));
    signal_rx(j,:) = filtered_signal_tx + noise(j,:);
    noise_uncoded(j,:) = sqrt(NoisePower_uncoded(j)/2).*(randn(1,length(signal_tx_uncoded)*M+N-1)+1i*randn(1,length(signal_tx_uncoded)*M+N-1));
    signal_rx_uncoded(j,:) = filtered_signal_tx_uncoded + noise_uncoded(j,:);
end


%%
% RRC Nyquist Filter RX
%-------------------------

filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*M+2*(N-1));
cropped_filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*M);
filtered_signal_rx_uncoded = zeros(length(EbN0),length(signal_tx_uncoded)*M+2*(N-1));
cropped_filtered_signal_rx_uncoded = zeros(length(EbN0),length(signal_tx_uncoded)*M);
for i =1:length(EbN0)
    filtered_signal_rx(i,:) = conv(signal_rx(i,:),fliplr(h_RRC));
    cropped_filtered_signal_rx(i,:) = filtered_signal_rx(i,N:end-(N-1));
    filtered_signal_rx_uncoded(i,:) = conv(signal_rx_uncoded(i,:),fliplr(h_RRC));
    cropped_filtered_signal_rx_uncoded(i,:) = filtered_signal_rx_uncoded(i,N:end-(N-1));
end

%%
% Downsampling
%-------------

downsampled_signal = zeros(length(EbN0),length(signal_tx));
downsampled_signal_uncoded = zeros(length(EbN0),length(signal_tx_uncoded));
for j = 1:length(EbN0)
    for i = 1:length(signal_tx)
        downsampled_signal(j,i)=cropped_filtered_signal_rx(j,1+M*(i-1));
    end
    for i = 1:length(signal_tx_uncoded)
        downsampled_signal_uncoded(j,i)=cropped_filtered_signal_rx_uncoded(j,1+M*(i-1));
    end
end

%%
%Demapping
%-----------

bits_rx = zeros(length(EbN0),length(bits_tx));
codedbits_rx = zeros(length(EbN0),length(codedbits_tx));
for i = 1:length(EbN0)
    if Nbps>1
        codedbits_rx(i,:) = demapping(downsampled_signal(i,:).',Nbps,"qam");
        bits_rx(i,:)=demapping(downsampled_signal_uncoded(i,:).',Nbps,"qam");
    else
        codedbits_rx(i,:) = demapping(real(downsampled_signal(i,:).'),Nbps,"pam");
        bits_rx(i,:)=demapping(real(downsampled_signal_uncoded(i,:).'),Nbps,"pam");
    end
end

%%
% Hard & Soft Decoding
%----------------

bits_rx_HD=zeros(length(EbN0),Nb);
bits_rx_SD=zeros(length(EbN0),Nb);
for i = 1:length(EbN0)
    for j = 1:BlockNb
        codeword = codedbits_rx(i,(j-1)*BlockSize/CodeRate+1:j*BlockSize/CodeRate);
        codeword2 = downsampled_signal(i,(j-1)*BlockSize/(CodeRate*Nbps)+1:j*BlockSize/(CodeRate*Nbps));
        correctedCodeword_HD=hardDecoding(codeword,H,10);
        correctedCodeword_SD=softDecoding(codeword2,H,N0(i)/2,10);
        bits_rx_HD(i,(j-1)*BlockSize+1:j*BlockSize)=correctedCodeword_HD(BlockSize+1:BlockSize/CodeRate);
        bits_rx_SD(i,(j-1)*BlockSize+1:j*BlockSize)=correctedCodeword_SD(BlockSize+1:BlockSize/CodeRate);
    end
end


%%
% BER
%----------

BER =zeros(1,length(EbN0));
BER_HD =zeros(1,length(EbN0));
BER_SD = zeros(1,length(EbN0));
for j = 1:length(EbN0)
    for i=1:Nb
        if(bits_rx(j,i) ~= bits_tx(1,i))
            BER(j) = BER(j)+1;
        end
        if(bits_rx_HD(j,i) ~= bits_tx(1,i))
            BER_HD(j) = BER_HD(j)+1;
        end
                if(bits_rx_SD(j,i) ~= bits_tx(1,i))
            BER_SD(j) = BER_SD(j)+1;
        end

    end
BER(j) = BER(j)/Nb;
BER_HD(j) = BER_HD(j)/Nb;
BER_SD(j) = BER_SD(j)/Nb;
end

AverageBER=AverageBER+BER;
AverageBER_HD=AverageBER_HD+BER_HD;
AverageBER_SD=AverageBER_SD+BER_SD;
end
AverageBER=AverageBER/AverageNb;
AverageBER_HD=AverageBER_HD/AverageNb;
AverageBER_SD=AverageBER_SD/AverageNb;

figure;
semilogy(EbN0,AverageBER)
hold on;
semilogy(EbN0,AverageBER_HD)
hold on;
semilogy(EbN0,AverageBER_SD)
hold off;
grid on;
title("BPSK (Nbps=1)");
legend('Uncoded','Hard Decoding (10 iter)','Soft Decoding (10 iter)');
xlabel("Eb/N0 [dB]");
ylabel("BER");