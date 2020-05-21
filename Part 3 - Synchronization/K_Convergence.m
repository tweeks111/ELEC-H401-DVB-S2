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
Nbps=4;                                        % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6) -> vector to compare 
CutoffFreq= 1e6;                                % CutOff Frequency of the Nyquist Filter
RollOff= 0.3;                                   % Roll-Off Factor
M= 50;                                          % Upsampling Factor
N = 16*M+1;                                     % Number of taps (ODD ONLY)
EbN0 = 1000;                                    % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
Tsymb= 1/(2*CutoffFreq);                        % Symbol Period
SymRate= 1/Tsymb;                               % Symbol Rate
Fs = SymRate*M;                                 % Sampling Frequency
BlockSize = 128;
BlockNb=30;
CodeRate = 1/2;
Nb= BlockSize*BlockNb;                          % Number of bits
timeShift = 20;
Fc = 2e9;
ppm = 0;
CFO = ppm*Fc*1e-6;
phase_offset_deg = 0;
phase_offset= phase_offset_deg*pi/180;
K=[0.01 0.02 0.05];
AverageNb= 50;
AverageTimeError = zeros(AverageNb,Nb/Nbps,length(K));


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
    signal_rx_sync_errors = signal_rx.*exp(1j*(2*pi*CFO.*t1+phase_offset));

    %%
    % RRC Nyquist Filter RX
    %-------------------------

    t2=((0:length(signal_tx)*M-1))*1/Fs;
    filtered_signal_rx = conv(signal_rx_sync_errors,h_RRC);

    %%
    % Time Shift
    %-----------------------
    cropped_filtered_signal_rx = filtered_signal_rx(N:end-(N-1));
    shifted_signal_rx=circshift(cropped_filtered_signal_rx,timeShift);
    
    %cropped_filtered_signal_rx = cropped_filtered_signal_rx.*exp(-1j*(2*pi*CFO.*t2));
                                                                  
    
    
    %%
    % Gardner
    %-------------
    downsampling_ratio=M/2;
    partial_downsampled_signal_rx = downsample(shifted_signal_rx,downsampling_ratio);
    downsampled_signal_rx_corrected = zeros(length(signal_tx),length(K));
    time_error = zeros(length(signal_tx),length(K));
    for k = 1:length(K)
        [downsampled_signal_rx_corrected(:,k),time_error(:,k)]=gardner(partial_downsampled_signal_rx,K(k),M/downsampling_ratio);
        AverageTimeError(avr,:,k)=time_error(:,k);
    end

end
    MeanTimeError = zeros(Nb/Nbps,length(K));
    VarianceTimeError = zeros(Nb/Nbps,length(K));
    for k = 1:length(K)
        MeanTimeError(:,k) = mean(AverageTimeError(:,:,k));
        MeanTimeError(:,k)=(MeanTimeError(:,k)+timeShift/M);
        VarianceTimeError(:,k) = std(AverageTimeError(:,:,k));
    end
colorVector = ['r','b','g','y','m','c','k',];
    Legend=cell(length(K));
for i = 1:length(K)
   vector = 1:25:Nb/Nbps;

   p(i) = plot(vector,MeanTimeError(vector,i),[colorVector(i) 'o-']);
   hold on;
   plot(vector,MeanTimeError(vector,i)-VarianceTimeError(vector,i),[colorVector(i) '--']);
   hold on;
   plot(vector,MeanTimeError(vector,i)+VarianceTimeError(vector,i),[colorVector(i) '--']);
   Legend{i}=['\kappa=' num2str(K(i))];
end

legend(p(1:length(K)),Legend(1:length(K)));

grid on;
xlabel("Symbols");
ylabel("Time error (mean±deviation)");
if(Nbps==1) 
    text='BPSK ';
elseif(Nbps==2) 
    text='QPSK ';
elseif(Nbps==4) 
    text='16QAM ';
else
    text ='64QAM ';
end

txt = {['M= ' num2str(M) ' | N= ' num2str(N)],['\beta= ' num2str(RollOff)],['t_0=' num2str(timeShift/M) 'T_{symb}'],['f_{symb}= ' num2str(SymRate*1e-6) 'MBd'],['CFO=' num2str(ppm) 'ppm | \phi_0=' num2str(phase_offset_deg) '°']};
annotation('textbox',[0.6,0.4,0.25,0.25],'String',txt,'BackgroundColor','white');

title([text,'(Nbps=',num2str(Nbps),')']);