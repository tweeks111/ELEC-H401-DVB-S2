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
Nbps= 2;                                        % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6) -> vector to compare 
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
ppm = 2;
CFO = ppm*Fc*1e-6;                              % Carrier Frequency Offset
phase_offset = 0;
time_shift=0;
K=[0 0.01 0.05];
AverageNb= 100;
AverageBER=zeros(1,length(EbN0));
AverageBERgardner=zeros(length(K),length(EbN0));


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
signal_rx_sync_errors=zeros(length(EbN0),size(signal_rx,2),length(CFO));

for i = 1:length(EbN0)
    signal_rx_sync_errors(i,:) = signal_rx(i,:).*exp(1j*(2*pi*CFO*t1+phase_offset));
end


%%
% RRC Nyquist Filter RX
%-------------------------

filtered_signal_rx = zeros(1,length(signal_tx)*M+2*(N-1));
cropped_filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*M);
t2=((0:length(signal_tx)*M-1))*1/Fs;

for i =1:length(EbN0)
    filtered_signal_rx = conv(signal_rx_sync_errors(i,:),fliplr(h_RRC));
    cropped_filtered_signal_rx(i,:) = filtered_signal_rx(N:end-(N-1));
end


%%
% Gardner
%---------------

signal_rx_corrected = zeros(length(EbN0),length(signal_tx),length(K));
time_error = zeros(length(EbN0),length(signal_tx),length(K));
for i =1:length(EbN0)
    for j = 1:length(K)
     [signal_rx_corrected(i,:,j),time_error(i,:,j)]=gardner(cropped_filtered_signal_rx(i,:),K(j),M);
    end
end


%%
% Downsampling
%-------------

downsampled_signal = zeros(length(EbN0),length(signal_tx));
for j = 1:length(EbN0)
    for i = 1:length(signal_tx)
        downsampled_signal(j,i)=cropped_filtered_signal_rx(j,1+M*(i-1));
    end
    downsampled_signal(j,:)=downsampled_signal(j,:);
end


%%
%Demapping
%-----------

bits_rx = zeros(length(EbN0),length(bits_tx));
bits_rx_gardner = zeros(length(EbN0),length(bits_tx),length(K));

for i = 1:length(EbN0)
    if Nbps>1
        bits_rx(i,:) = demapping(downsampled_signal(i,:).',Nbps,"qam");
        for j = 1:length(K)
            bits_rx_gardner(i,:,j) = demapping(signal_rx_corrected(i,:,j).',Nbps,"qam");
            
        end
    else
        bits_rx(i,:) = demapping(real(downsampled_signal(i,:).'),Nbps,"pam");
        for j = 1:length(K)
            bits_rx_gardner(i,:,j) = demapping(real(signal_rx_corrected(i,:,j).'),Nbps,"pam");
        end
    end
end

%%
% BER
%----------

BER =zeros(1,length(EbN0));
BER_gardner = zeros(length(K),length(EbN0));
for j = 1:length(EbN0)
    for i=1:Nb
        if(bits_rx(j,i) ~= bits_tx(1,i))
            BER(j) = BER(j)+1/Nb;
        end
        for k = length(K)
            if(bits_rx_gardner(j,i,k)~=bits_tx(1,i))
                BER_gardner(k,j)=BER_gardner(k,j)+1;
            end
        end
    end
end


AverageBER=AverageBER+BER;
AverageBERgardner=AverageBERgardner+BER_gardner;
end


AverageBER=AverageBER/AverageNb;


figure;
for k = 1:length(K)
        label=['K=' num2str(K(k))]; 
        semilogy(EbN0,AverageBERgardner(k,:),'DisplayName',label);
end
hold on;
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