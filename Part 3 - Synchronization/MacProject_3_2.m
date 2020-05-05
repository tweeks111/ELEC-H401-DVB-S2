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
addpath('../Part 2 - LDPC');
%------Parameters------%
Nbps= 4;                                        % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6) -> vector to compare 
CutoffFreq= 1e6;                                % CutOff Frequency of the Nyquist Filter
RollOff= 0.3;                                   % Roll-Off Factor
M= 16;                                          % Upsampling Factor
N = 101;                                        % Number of taps (ODD ONLY)
EbN0 = -2:1:20;                                 % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
Tsymb= 1/(2*CutoffFreq);                        % Symbol Period
SymRate= 1/Tsymb;                               % Symbol Rate
Fs = SymRate*M;                                 % Sampling Frequency
BlockSize = 128;
BlockNb=6;
CodeRate = 1/2;
Nb= BlockSize*BlockNb;                          % Number of bits
AverageNb=50;
AverageBER=zeros(4,length(EbN0));
H0 = makeLdpc(BlockSize, BlockSize/CodeRate,0,1,3);
Fc = 2e9;
ppm = 40;
CFO = ppm*Fc*1e-6;                              % Carrier Frequency Offset
phase_offset = 0;
time_shift=0;


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
else
        signal_tx = mapping(codedbits_tx,Nbps,'pam').';         % Symbols sequence at transmitter   
end

%%
% Upsampling
%-----------------

upsampled_signal = zeros(1,length(signal_tx)*M);
for i = 1:length(signal_tx)
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
Eb = SignalEnergy/(2*Nb/CodeRate);

N0 = Eb./(10.^(EbN0/10));
NoisePower = 2*N0*Fs;

noise = zeros(length(EbN0),length(signal_tx)*M+N-1);
signal_rx = zeros(length(EbN0),length(signal_tx)*M+N-1);

for j = 1:length(EbN0)
    noise(j,:) = sqrt(NoisePower(j)/2).*(randn(1,length(signal_tx)*M+N-1)+1i*randn(1,length(signal_tx)*M+N-1));
    signal_rx(j,:) = filtered_signal_tx + noise(j,:);
end

%%
% CFO & Carrier Phase Error
%--------------------

t = ((0:length(signal_rx)-1)-(N-1)/2)*1/Fs;
signal_rx_sync_errors=zeros(length(EbN0),length(signal_rx));
for i = 1:length(EbN0)
    signal_rx_sync_errors(i,:) = signal_rx(i,:).*exp(1j*(2*pi*CFO*t+phase_offset));
end

%%
% RRC Nyquist Filter RX
%-------------------------

filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*M+2*(N-1));
cropped_filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*M);
filtered_signal_rx_sync_errors = zeros(length(EbN0),length(signal_tx)*M+2*(N-1));
cropped_filtered_signal_rx_sync_errors = zeros(length(EbN0),length(signal_tx)*M);
t=((0:length(cropped_filtered_signal_rx)-1))*1/Fs;
for i =1:length(EbN0)
    filtered_signal_rx(i,:) = conv(signal_rx(i,:),fliplr(h_RRC));
    cropped_filtered_signal_rx(i,:) = filtered_signal_rx(i,N:end-(N-1));
    filtered_signal_rx_sync_errors(i,:) = conv(signal_rx_sync_errors(i,:),fliplr(h_RRC));
    cropped_filtered_signal_rx_sync_errors(i,:) = filtered_signal_rx_sync_errors(i,N:end-(N-1));
    cropped_filtered_signal_rx_sync_errors(i,:)=cropped_filtered_signal_rx_sync_errors(i,:).*exp(-1j*(2*pi*CFO*t+phase_offset));
end
 

%%
% Downsampling
%-------------

downsampled_signal = zeros(length(EbN0),length(signal_tx));
downsampled_signal_sync_errors = zeros(length(EbN0),length(signal_tx));
for j = 1:length(EbN0)
    for i = 1:length(signal_tx)
        downsampled_signal(j,i)=cropped_filtered_signal_rx(j,1+M*(i-1));
        downsampled_signal_sync_errors(j,i)=cropped_filtered_signal_rx_sync_errors(j,1+M*(i-1));
    end
end



%%
%Demapping
%-----------

codedbits_rx = zeros(length(EbN0),length(codedbits_tx));
codedbits_rx_sync_errors = zeros(length(EbN0),length(codedbits_tx));
for i = 1:length(EbN0)
    if Nbps>1
        codedbits_rx(i,:) = demapping(downsampled_signal(i,:).',Nbps,"qam");
        codedbits_rx_sync_errors(i,:) = demapping(downsampled_signal_sync_errors(i,:).',Nbps,"qam");
    else
        codedbits_rx(i,:) = demapping(real(downsampled_signal(i,:).'),Nbps,"pam");
        codedbits_rx_sync_errors(i,:) = demapping(real(downsampled_signal_sync_errors(i,:).'),Nbps,"pam");
    end
end

%%
% Hard & Soft Decoding
%----------------

bits_rx_uncoded=zeros(length(EbN0),Nb);
bits_rx_uncoded_sync_errors=zeros(length(EbN0),Nb);
% bits_rx_HD=zeros(length(EbN0),Nb);
% bits_rx_SD=zeros(length(EbN0),Nb);
for i = 1:length(EbN0)
    for j = 1:BlockNb
        codeword = codedbits_rx(i,(j-1)*BlockSize/CodeRate+1:j*BlockSize/CodeRate);
        codeword_sync_errors = codedbits_rx_sync_errors(i,(j-1)*BlockSize/CodeRate+1:j*BlockSize/CodeRate);
        
%        coded_symbols = downsampled_signal(i,(j-1)*BlockSize/(CodeRate*Nbps)+1:j*BlockSize/(CodeRate*Nbps));
%        correctedCodeword_SD=softDecoding(coded_symbols,H,N0(i)/2,10);
%        correctedCodeword_HD=hardDecoding(codeword,H,10);
%        bits_rx_HD(i,(j-1)*BlockSize+1:j*BlockSize)=correctedCodeword_HD(BlockSize+1:BlockSize/CodeRate);
%        bits_rx_SD(i,(j-1)*BlockSize+1:j*BlockSize)=correctedCodeword_SD(BlockSize+1:BlockSize/CodeRate);
        bits_rx_uncoded(i,(j-1)*BlockSize+1:j*BlockSize)=codeword(BlockSize+1:BlockSize/CodeRate);
        bits_rx_uncoded_sync_errors(i,(j-1)*BlockSize+1:j*BlockSize)=codeword_sync_errors(BlockSize+1:BlockSize/CodeRate);
    end
end


%%
% BER
%----------

BER =zeros(4,length(EbN0));
for j = 1:length(EbN0)
    for i=1:Nb
        if(bits_rx_uncoded(j,i) ~= bits_tx(1,i))
            BER(1,j) = BER(1,j)+1/Nb;
        end
        if(bits_rx_uncoded_sync_errors(j,i)~= bits_tx(1,i))     
            BER(2,j) = BER(2,j)+1/Nb;
        end
%         if(bits_rx_HD(j,i)~= bits_tx(1,i))
%             BER(3,j) = BER(3,j)+1/Nb;
%         end
%         if(bits_rx_SD(j,i)~= bits_tx(1,i))
%             BER(4,j) = BER(4,j)+1/Nb;
%         end
    end
end

AverageBER=AverageBER+BER;
end


AverageBER=AverageBER/AverageNb;

figure;
semilogy(EbN0,AverageBER(1,:))
hold on;
semilogy(EbN0,AverageBER(2,:))
hold on;
% semilogy(EbN0,AverageBER(3,:))
% hold on;
% semilogy(EbN0,AverageBER(4,:))
% hold off;
grid on;
title("16QAM (Nbps=4)");
legend('No CFO','CFO - ppm=10');
xlabel("Eb/N0 [dB]");
ylabel("BER");