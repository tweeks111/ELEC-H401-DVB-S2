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
M= 50;                                           % Upsampling Factor
N = 16*M+1;                                            % Number of taps (ODD ONLY)
EbN0 = 1000;                                 % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
Tsymb= 1/(2*CutoffFreq);                        % Symbol Period
SymRate= 1/Tsymb;                               % Symbol Rate
Fs = SymRate*M;                                 % Sampling Frequency
BlockSize = 128;
BlockNb=30;
CodeRate = 1/2;
Nb= BlockSize*BlockNb;                          % Number of bits
%H0 = makeLdpc(BlockSize, BlockSize/CodeRate,0,1,3);
Fc = 2e9;
ppm = [0 2 50 100 150 175 200];
CFO = ppm*Fc*1e-6;                              % Carrier Frequency Offset
phase_offset_deg = 0;
phase_offset= phase_offset_deg*pi/180;
K=0.05;
timeShift=20;
AverageNb= 50;
AverageTimeError = zeros(AverageNb,Nb/Nbps,length(CFO));

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

    noise = sqrt(NoisePower/2).*(randn(1,length(signal_tx)*M+N-1)+1i*randn(1,length(signal_tx)*M+N-1));
    signal_rx = filtered_signal_tx + noise;

    %%
    % CFO & Carrier Phase Error
    %--------------------

    t1 = ((0:size(signal_rx,2)-1))*1/Fs;
    signal_rx_sync_errors=zeros(size(signal_rx,2),length(CFO));
    for k=1:length(CFO)
       signal_rx_sync_errors(:,k) = signal_rx.*exp(1j*(2*pi*CFO(k).*t1+phase_offset));
    end

    %%
    % RRC Nyquist Filter RX
    %-------------------------

    filtered_signal_rx = zeros(1,length(signal_tx)*M+2*(N-1));
    cropped_filtered_signal_rx = zeros(length(signal_tx)*M,length(CFO));
    t2=((0:length(signal_tx)*M-1))*1/Fs;
    for k = 1:length(CFO)
            filtered_signal_rx = conv(signal_rx_sync_errors(:,k),fliplr(h_RRC));
            cropped_filtered_signal_rx(:,k) = filtered_signal_rx(N:end-(N-1));
            %cropped_filtered_signal_rx(:,k) = cropped_filtered_signal_rx(i,:,k).*exp(-1j*2*pi*CFO(k).*t2);
                                                                  %           /\          %
    end                                                                          %  To observe ISI only  %
    
    %%
    % Time Shift
    %-----------------------
    downsampling_ratio=M/2;
    shifted_signal_rx = zeros(length(signal_tx)*M,length(CFO));
    partial_downsampled_signal_rx = zeros(length(signal_tx)*M/downsampling_ratio,length(CFO));
    
    for k =1:length(CFO)
            shifted_signal_rx(:,k)=circshift(cropped_filtered_signal_rx(:,k),timeShift);
            partial_downsampled_signal_rx(:,k) = downsample(shifted_signal_rx(:,k),downsampling_ratio);
    end                                                          
    
    
    %%
    % Gardner
    %-------------

    downsampled_signal_rx_corrected = zeros(length(signal_tx),length(CFO));
    time_error = zeros(length(signal_tx),length(CFO));
    for k = 1:length(CFO)
        [downsampled_signal_rx_corrected(:,k),time_error(:,k)]=gardner(partial_downsampled_signal_rx(:,k),K,M/downsampling_ratio);
        AverageTimeError(avr,:,k)=time_error(:,k);
    end
    
end

MeanTimeError = zeros(Nb/Nbps,length(CFO));
    for k = 1:length(CFO)
        MeanTimeError(:,k) = mean(AverageTimeError(:,:,k));
        MeanTimeError(:,k)=(MeanTimeError(:,k)+timeShift/M);
    end
colorVector = ['r','b','g','m','c','k','y'];
    Legend=cell(length(CFO));
    
for i = 1:length(CFO)
   vector = 1:25:Nb/Nbps;

   p(i) = plot(vector,MeanTimeError(vector,i),[colorVector(i) 'o-']);
   hold on;
   Legend{i}=['CFO=' num2str(ppm(i)) 'ppm'];
end

legend(p(1:length(CFO)),Legend(1:length(CFO)));

grid on;

legend('show');
xlabel("Symbols");
ylabel("Time error (mean)");
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
% txt = {['#taps= ' num2str(N)],['RollOff= ' num2str(RollOff)],['M= ' num2str(M)],['SymRate= ' num2str(SymRate*1e-6) 'MBd'],['Phase Offset=' num2str(phase_offset_deg) '°']};
% annotation('textbox',[0.2,0.2,0.22,0.22],'String',txt,'BackgroundColor','white');

title([text,'(Nbps=',num2str(Nbps),')']);