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
M= 8;                                          % Upsampling Factor
N = 51;                                        % Number of taps (ODD ONLY)
EbN0 = -5:1:20;                                 % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
Tsymb= 1/(2*CutoffFreq);                        % Symbol Period
SymRate= 1/Tsymb;                               % Symbol Rate
Fs = SymRate*M;                                 % Sampling Frequency
BlockSize = 128;
BlockNb=6;
CodeRate = 1/2;
Nb= BlockSize*BlockNb;                          % Number of bits
%H0 = makeLdpc(BlockSize, BlockSize/CodeRate,0,1,3);
Fc = 2e9;
ppm = 0;
CFO = ppm*Fc*1e-6;                              % Carrier Frequency Offset
phase_offset_deg = [1 2 3 5 10 30];
phase_offset = phase_offset_deg*pi/180;
AverageNb= 100;
AverageBER=zeros(length(phase_offset),length(EbN0));


for avr = 1:AverageNb
disp(avr);
%%
% Bit Generation
%------------------------

bits_tx = randi(2,1,Nb)-1;               % bits_tx = Binary sequence


%%
% Mapping
%------------------------

if Nbps>1
        signal_tx = mapping(bits_tx.',Nbps,'qam').';         % Symbols sequence at transmitter
else
        signal_tx = mapping(bits_tx.',Nbps,'pam').';         % Symbols sequence at transmitter   
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
Eb = SignalEnergy/(2*Nb);

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

t1 = ((0:size(signal_rx,2)-1))*1/Fs;
signal_rx_sync_errors=zeros(length(EbN0),size(signal_rx,2),length(phase_offset));
for k=1:length(phase_offset)
    for i = 1:length(EbN0)
        signal_rx_sync_errors(i,:,k) = signal_rx(i,:).*exp(1j*(2*pi*CFO*t1+phase_offset(k)));
    end
end

%%
% RRC Nyquist Filter RX
%-------------------------

filtered_signal_rx = zeros(1,length(signal_tx)*M+2*(N-1));
cropped_filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*M,length(phase_offset));
t2=((0:length(signal_tx)*M-1))*1/Fs;
for k = 1:length(phase_offset)
    for i =1:length(EbN0)
        filtered_signal_rx = conv(signal_rx_sync_errors(i,:,k),fliplr(h_RRC));
        cropped_filtered_signal_rx(i,:,k) = filtered_signal_rx(N:end-(N-1));
        cropped_filtered_signal_rx(i,:,k) = cropped_filtered_signal_rx(i,:,k).*exp(-1j*(2*pi*CFO.*t2));
    end                                                                      %           /\
end                                                                          %  To observe ISI only

%%
% Downsampling
%-------------

downsampled_signal = zeros(length(EbN0),length(signal_tx),length(phase_offset));
for k = 1:length(phase_offset)
    for j = 1:length(EbN0)
        for i = 1:length(signal_tx)
            downsampled_signal(j,i,k)=cropped_filtered_signal_rx(j,1+M*(i-1),k);
        end
        downsampled_signal(j,:,k)=downsampled_signal(j,:,k);
    end
end

%%
%Demapping
%-----------

bits_rx = zeros(length(EbN0),length(bits_tx),length(phase_offset));
for k = 1:length(phase_offset)
    for i = 1:length(EbN0)
        if Nbps>1
            bits_rx(i,:,k) = demapping(downsampled_signal(i,:,k).',Nbps,"qam");
        else
            bits_rx(i,:,k) = demapping(real(downsampled_signal(i,:,k).'),Nbps,"pam");
        end
    end
end
%%
% BER
%----------

BER =zeros(length(phase_offset),length(EbN0));
for k = 1:length(phase_offset)
    for j = 1:length(EbN0)
        for i=1:Nb
            if(bits_rx(j,i,k) ~= bits_tx(1,i))
                BER(k,j) = BER(k,j)+1/Nb;
            end
        end
    end
end

AverageBER=AverageBER+BER;
end


AverageBER=AverageBER/AverageNb;


figure;
for k = 1:length(phase_offset_deg)
    if(phase_offset_deg(k)==0)
        label='No Phase Offset';
    else
        label=['Phase Offset =' num2str(phase_offset_deg(k)) '°']; 
    end
semilogy(EbN0,AverageBER(k,:),'DisplayName',label);
hold on;
end
grid on;

legend('show');
xlabel("Eb/N0 [dB]");
ylabel("BER");
if(Nbps==1) 
    text='BPSK ';
elseif(Nbps==2) 
    text='QPSK ';
elseif(Nbps==4) 
    text='16QAM ';
else
    text ='64QAM ';
end

txt = {['#taps= ' num2str(N)],['RollOff= ' num2str(RollOff)],['M= ' num2str(M)],['SymRate: ' num2str(SymRate*1e-6) 'MBd']};
annotation('textbox',[0.2,0.2,0.22,0.18],'String',txt,'BackgroundColor','white');

title([text,'(Nbps=',num2str(Nbps),')']);