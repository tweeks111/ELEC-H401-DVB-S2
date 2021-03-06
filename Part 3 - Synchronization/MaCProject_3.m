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
Nbps= 1;                                        % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6) -> vector to compare 
CutoffFreq= 1e6;                                % CutOff Frequency of the Nyquist Filter
RollOff= 0.3;                                   % Roll-Off Factor
M= 50;                                           % Upsampling Factor
N = 16*M+1;                                            % Number of taps (ODD ONLY)
EbN0 = 10;                                 % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
Tsymb= 1/(2*CutoffFreq);                        % Symbol Period
SymRate= 1/Tsymb;                               % Symbol Rate
Fs = SymRate*M;                                 % Sampling Frequency
BlockSize = 128;
BlockNb=10;
CodeRate = 1/2;
Nb= BlockSize*BlockNb;                          % Number of bits
%H0 = makeLdpc(BlockSize, BlockSize/CodeRa te,0,1,3);
Fc = 2e9;
ppm = 10;
CFO = ppm*Fc*1e-6;                              % Carrier Frequency Offset
phase_offset_deg = 0;
phase_offset= phase_offset_deg*pi/180;
K=0.05;
timeShift=10;
AverageNb= 50;
AverageTimeError = zeros(AverageNb,Nb/Nbps);
AverageBER=0;

pilot_pos = 34;
pilot_size = 20;
avgWindow_size = 8;

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

    unuseful = signal_tx(1:pilot_pos-1);
    pilot = signal_tx(pilot_pos : pilot_pos+pilot_size-1);
    symbols = signal_tx(pilot_pos+pilot_size : end);

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

    filtered_signal_rx = conv(signal_rx_sync_errors,fliplr(h_RRC));
    cropped_filtered_signal_rx = filtered_signal_rx(N:end-(N-1));
    
    %%
    % Time Shift
    %-----------------------
    
    downsampling_ratio=M/2;

    shifted_signal_rx=circshift(cropped_filtered_signal_rx,timeShift);
    partial_downsampled_signal_rx = downsample(shifted_signal_rx,downsampling_ratio);                                                  
    
    %%
    % Gardner
    %-------------

    [downsampled_signal_rx_corrected,time_error]=gardner(partial_downsampled_signal_rx,K,M/downsampling_ratio);
    AverageTimeError(avr,:)=time_error;
    
    %%
    % Data acquisition
    %-----------------------------------

    [toa, est_CFO] = dataAcquisition(downsampled_signal_rx_corrected,pilot,avgWindow_size, Tsymb);
    disp(["ToA :" num2str(toa)]);
    disp(["Estimated CFO : " num2str(est_CFO)]);
    t2=((0:length(signal_tx)-1))*Tsymb;
    cfo_corrected_signal_rx = downsampled_signal_rx_corrected.*exp(-1j*2*pi*est_CFO*t2);
    
    %%
    %Demapping
    %-----------
    
    if Nbps>1
        bits_rx = demapping(cfo_corrected_signal_rx.',Nbps,"qam");    
    else
        bits_rx = demapping(real(cfo_corrected_signal_rx.'),Nbps,"pam");
    end
    
    %%
    % BER
    %----------

    BER =0;
    
    for i=1:Nb
        if(bits_rx(i) ~= bits_tx(i))
            BER = BER+1/Nb;
        end
    end
    
  
    AverageBER=AverageBER+BER/AverageNb;
end
    disp(["BER = " AverageBER]);

%%
% PLOTS
%------------
if(Nbps==1) 
    text='BPSK ';
elseif(Nbps==2) 
    text='QPSK ';
elseif(Nbps==4) 
    text='16QAM ';
else
    text ='64QAM ';
end

% Plot error convergence
%------------------------
MeanTimeError = mean(AverageTimeError);
MeanTimeError=(MeanTimeError+timeShift/M);

vector = 1:25:Nb/Nbps;
plot(vector,MeanTimeError(vector),'ro-');
hold on;
legend(['CFO=' num2str(ppm) 'ppm']);

grid on;
legend('show');
xlabel("Symbols");
ylabel("Time error (mean)");
% 
% txt = {['#taps= ' num2str(N)],['RollOff= ' num2str(RollOff)],['M= ' num2str(M)],['SymRate= ' num2str(SymRate*1e-6) 'MBd'],['Phase Offset=' num2str(phase_offset_deg) '°']};
% annotation('textbox',[0.2,0.2,0.22,0.22],'String',txt,'BackgroundColor','white');

title([text,'(Nbps=',num2str(Nbps),')']);

