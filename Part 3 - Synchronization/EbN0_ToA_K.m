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
Nbps=2;                                        % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6) -> vector to compare 
CutoffFreq= 1e6;                                % CutOff Frequency of the Nyquist Filter
RollOff= 0.3;                                   % Roll-Off Factor
M= 50;                                          % Upsampling Factor
N = 16*M+1;                                     % Number of taps (ODD ONLY)
EbN0 = 0:2:16;                                    % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
Tsymb= 1/(2*CutoffFreq);                        % Symbol Period
SymRate= 1/Tsymb;                               % Symbol Rate
Fs = SymRate*M;                                 % Sampling Frequency
BlockSize = 128;
BlockNb=10;
CodeRate = 1/2;
Nb= BlockSize*BlockNb;                          % Number of bits
timeShift = 0;
Fc = 2e9;
ppm = 0;
CFO = ppm*Fc*1e-6;
phase_offset_deg = 0;
phase_offset= phase_offset_deg*pi/180;
K=0.05;
AverageNb= 100;


ToA = 34;
pilot_size = 20;
avgWindow_size = [1 8 16];
AverageTimeError = zeros(AverageNb,length(EbN0),length(avgWindow_size));

for avr = 1:AverageNb
    %%
    % Bit Generation
    %------------------------
    
    disp(avr);
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
    %Divide msg into unuseful data, pilot, symbols
    %----------------------------------------------

    unuseful = signal_tx(1:ToA-1);
   
    pilot = signal_tx(ToA : ToA+pilot_size-1);
    symbols = signal_tx(ToA+pilot_size: end);
    
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
    % RRC Nyquist Filter RX
    %-------------------------

    filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*M+2*(N-1));
    for i =1:length(EbN0)
        filtered_signal_rx(i,:) = conv(signal_rx(i,:),fliplr(h_RRC));
    end                                                                      

    %%
    % Time Shift
    %-----------------------
    shifted_signal_rx = zeros(length(EbN0),length(filtered_signal_rx));
    cropped_filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*M);
    t2=((0:length(signal_tx)*M-1))*1/Fs;
    for i = 1:length(EbN0)
        shifted_signal_rx(i,:)=circshift(filtered_signal_rx(i,:),timeShift);
        cropped_filtered_signal_rx(i,:) = shifted_signal_rx(i,N:end-(N-1));
        %cropped_filtered_signal_rx(i,:,k) = cropped_filtered_signal_rx(i,:,k).*exp(-1j*(2*pi*CFO.*t2));
    end

%     %%
%     % Gardner
%     %-------------
%     downsampling_ratio=M/2;
%     partial_downsampled_signal_rx = zeros(length(EbN0),size(cropped_filtered_signal_rx,2)/downsampling_ratio);
%     downsampled_signal_rx_corrected = zeros(length(EbN0),length(signal_tx));
%     time_error = zeros(length(EbN0),length(signal_tx));
%     for i = 1:length(EbN0)
%         partial_downsampled_signal_rx(i,:) = downsample(cropped_filtered_signal_rx(i,:),downsampling_ratio);
%         [downsampled_signal_rx_corrected(i,:),time_error(i,:)]=gardner(partial_downsampled_signal_rx(i,:),K,M/downsampling_ratio);
%     end
    
    %%
    % Downsample
    %-----------------
    
    downsampled_signal_rx = zeros(length(EbN0),length(signal_tx));
    for i = 1:length(EbN0)
        downsampled_signal_rx(i,:) = downsample(cropped_filtered_signal_rx(i,:),M);
    end
    
    %%
    % Data acquisition
    %-----------------------------------

    for i = 1:length(EbN0)
            [est_ToA, est_CFO] = dataAcquisition(downsampled_signal_rx(i,:),pilot,avgWindow_size(1), Tsymb);
            AverageTimeError(avr,i,1)=abs(est_ToA-ToA);
            [est_ToA, est_CFO] = dataAcquisition(downsampled_signal_rx(i,:),pilot,avgWindow_size(2), Tsymb);
            AverageTimeError(avr,i,2)=abs(est_ToA-ToA);
            [est_ToA, est_CFO] = dataAcquisition(downsampled_signal_rx(i,:),pilot,avgWindow_size(3), Tsymb);
            AverageTimeError(avr,i,3)=abs(est_ToA-ToA);
    end

    
end

VarianceTimeError = zeros(length(EbN0),length(avgWindow_size));
for i=1:length(avgWindow_size)
    for j = 1:length(EbN0)
        matrix = AverageTimeError(:,j,i).';
        VarianceTimeError(j,i) = std(matrix);
    end
end 

plot(EbN0,VarianceTimeError(:,1),'ro-');
hold on;
plot(EbN0,VarianceTimeError(:,2),'bo-');
hold on;
plot(EbN0,VarianceTimeError(:,3),'go-');
hold on;
grid on;
xlabel("EbN0");
ylabel("Time error stdev");%±deviation

if(Nbps==1) 
    text='BPSK ';
elseif(Nbps==2) 
    text='QPSK ';
elseif(Nbps==4) 
    text='16QAM ';
else
    text ='64QAM ';
end
% 
% txt = {['M= ' num2str(M) ' | N= ' num2str(N)],['\beta= ' num2str(RollOff)],['t_0=' num2str(timeShift/M) 'T_{symb}'],['f_{symb}= ' num2str(SymRate*1e-6) 'MBd'],['CFO=' num2str(ppm) 'ppm | \phi_0=' num2str(phase_offset_deg) '°']};
% annotation('textbox',[0.4,0.6,0.25,0.25],'String',txt,'BackgroundColor','white');

title([text,'(Nbps=',num2str(Nbps),')']);