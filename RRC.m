function [h_RRC,H_RRC] = RRC(Fs,Tsymb,N,RollOff,Nbps,AverageNb,USF)
% 
% 
% 
% 

df = (1/N)*Fs;                                 % Delta_f : 1 tap = df [Hz]
fmax = df*(N-1)/2;
fvector = linspace(-fmax,fmax,N);  

i=1;
H_RC=zeros(1,N);
for f = -fmax:df:fmax
    if (abs(f)<=(1-RollOff)/(2*Tsymb) )
       H_RC(i)=Tsymb;
    elseif(abs(f)<=(1+RollOff)/(2*Tsymb))
       H_RC(i)=Tsymb*(1+cos(pi*Tsymb*(abs(f)-(1-RollOff)/(2*Tsymb))/RollOff))/2;
    else
       H_RC(i)=0;
    end
    i=i+1;
end

H_RC = ifftshift(H_RC);
H_RRC = sqrt(H_RC);       
h_RC = ifft(H_RC,'symmetric');
h_RRC = ifft(H_RRC,'symmetric');
h_RRC = fftshift(h_RRC/sqrt(h_RC(1)));
h_RC = fftshift(h_RC/h_RC(1)); 
 
dt = 1/Fs;   
tvector = (-(N-1)/2-1:(N-1)/2-1)*dt;

if(length(Nbps)==1&&AverageNb==1)
    figure;
    plot(fvector,fftshift(H_RRC),'r-',fvector,fftshift(H_RRC),'b*');
    hold off;
    figure("Name","Impulse responses of the RC and RRC filters");
    p=zeros(1,4);
    p(1)=plot(tvector*10^6,h_RC,'r-');hold on;
    p(3)=plot(tvector*10^6,h_RRC,'g-');hold on
    p(5)=plot(tvector((N-1)/2+1:USF:N)*10^6,h_RC((N-1)/2+1:USF:N),'bo');hold on;
    p(6)=plot(tvector((N-1)/2+1:-USF:1)*10^6,h_RC((N-1)/2+1:-USF:1),'bo');
    hold off
    legend([p(1) p(3)],"h_{RC}","h_{RRC}");
    xlabel('Time [µs]')
end
end

