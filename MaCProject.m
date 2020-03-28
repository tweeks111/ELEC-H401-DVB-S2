%-------------------------------------%
%    Modulation and Coding Project    %
%-------------------------------------%
%   Authors : Theo LEPOUTTE           %
%             John ANTOUN             %
%                                     %
%   Date : March 16, 2020             %
%-------------------------------------%
clc;clear;close all;
%------Parameters------%
Nb= 6000;                  % Number of bits  
Nbps= [1 2 4 6];            % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6) -> vector to compare 
CutoffFreq= 1000000;        % CutOff Frequency of the Nyquist Filter
RollOff= 0.3;               % Roll-Off Factor
USF= 8;                     % Upsampling Factor
N = 101;                    % Number of taps (ODD ONLY)
EbN0 = 0:1:15;              % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
AverageNb = 50;             % Number of iteration to average the BER 
Tsymb= 1/(2*CutoffFreq);    % Symbol Period
SymRate= 1/Tsymb;           % Symbol Rate
Fs = USF*SymRate;           % Sampling Frequency
AverageBER=zeros(length(EbN0),length(Nbps));

% To display graphs, put AverageNb=1, Nbps = int and EbN0 = int

%=============================================%
% Bit Generation
%------------------------


for nb = 1:AverageNb

bits_tx = randi(2,1,Nb)-1;               % bits_tx = Binary sequence


for nbps = 1:length(Nbps)
    
% Mapping
%------------------------
if Nbps(nbps)>1
    signal_tx = mapping(bits_tx.',Nbps(nbps),'qam').';         % Symbols sequence at transmitter
else
    signal_tx = mapping(bits_tx.',Nbps(nbps),'pam').';         
end

% Upsampling
%-----------------------

upsampled_signal = zeros(1,Nb/Nbps(nbps)*USF);
for i = 1:Nb/Nbps(nbps)
    upsampled_signal(1+USF*(i-1))=signal_tx(i);
    for j = 2:USF
        upsampled_signal(j+USF*(i-1))=0;
    end
end

% RRC Nyquist Filter TX
%-------------------------

[h_RRC,H_RRC] =  RRC(Fs,Tsymb,N,RollOff,Nbps,AverageNb,USF);
filtered_signal_tx = conv(upsampled_signal,h_RRC);

if (length(EbN0)==1 && AverageNb==1)
    figure("Name","TX signal")
    subplot(1,2,1)
    plot(upsampled_signal,'ro')
    title("Upsampled TX signal")
    subplot(1,2,2)
    plot(filtered_signal_tx,"r.")
    title("Filtered TX signal")
end

% Noise
%-----------------

SignalEnergy = (trapz(abs(filtered_signal_tx).^2))*(1/Fs);
Eb = SignalEnergy/(2*Nb);

N0 = Eb./(10.^(EbN0/10));
NoisePower = 2*N0*Fs;

noise = zeros(length(EbN0),Nb/Nbps(nbps)*USF+N-1);
signal_rx = zeros(length(EbN0),Nb/Nbps(nbps)*USF+N-1);
for j = 1:length(EbN0)
    noise(j,:) = sqrt(NoisePower(j)/2).*(randn(1,Nb/Nbps(nbps)*USF+N-1)+1i*randn(1,Nb/Nbps(nbps)*USF+N-1));
    signal_rx(j,:) = filtered_signal_tx + noise(j,:);
end
% RRC Nyquist Filter RX
%-------------------------

filtered_signal_rx = zeros(length(EbN0),Nb/Nbps(nbps)*USF+2*(N-1));
cropped_filtered_signal_rx = zeros(length(EbN0),Nb/Nbps(nbps)*USF);
for i =1:length(EbN0)
    filtered_signal_rx(i,:) = conv(signal_rx(i,:),fliplr(h_RRC));
    cropped_filtered_signal_rx(i,:) = filtered_signal_rx(i,N:end-(N-1));
end

% Downsampling
%-------------

downsampled_signal = zeros(length(EbN0),Nb/Nbps(nbps));
for j = 1:length(EbN0)
    for i = 1:Nb/Nbps(nbps)
        downsampled_signal(j,i)=cropped_filtered_signal_rx(j,1+USF*(i-1));
    end
end

if (length(EbN0)==1 && AverageNb==1)
    figure("Name","RX signal");
    subplot(1,3,1);
    plot(signal_rx,"r.")
    title("Noised RX signal");
    subplot(1,3,2);
    plot(cropped_filtered_signal_rx,"r.")
    title("Filtered RX signal");
    subplot(1,3,3);
    plot(downsampled_signal,"r.");
    title("Cropped RX signal");
end
% Demapping
%-----------

bits_rx = zeros(length(EbN0),Nb);
for j = 1:length(EbN0)
    if Nbps(nbps)>1
        bits_rx(j,:) = demapping(downsampled_signal(j,:).',Nbps(nbps),"qam");
    else
        bits_rx(j,:) = demapping(real(downsampled_signal(j,:).'),Nbps(nbps),"pam");
    end
end

% BER
%-----

BER =zeros(length(EbN0),1);
for j = 1:length(EbN0)
    for i=1:Nb
        if(bits_rx(j,i) ~= bits_tx(i))
            BER(j,1) = BER(j,1)+1;
        end
    end
BER(j,1) = BER(j,1)/Nb;
end
AverageBER(:,nbps) = AverageBER(:,nbps) + BER(:,1)/AverageNb;
end
end

figure("Name","BER");
if(length(EbN0)>1)
    for nbps = 1:length(Nbps)
        if(Nbps(nbps)==1) 
            text=' 2-PAM';
        elseif(Nbps(nbps)==2) 
            text=' 4-QAM';
        elseif(Nbps(nbps)==4) 
            text='16-QAM';
        else
            text ='64-QAM';
        end
            
        
        semilogy(EbN0,AverageBER(:,nbps),'DisplayName', [text ' (Nbps=' num2str(Nbps(nbps)) ')']);
       
        hold on;
    end
end
hold off;
grid on;
legend('show');
xlabel("Eb/N0 [dB]");
ylabel("BER");

