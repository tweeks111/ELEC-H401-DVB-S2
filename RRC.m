function [h_RRC,H_RRC] = RRC(Fs,Tsymb,N,RollOff)
% 
% 
% 
% 

df = Fs/N;                                 % Delta_f : 1 tap = df [Hz]
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

H_RC = fftshift(H_RC);
H_RRC = sqrt(H_RC);       % normalizing RC filter bef
h_RC = ifft(H_RC,'symmetric');
h_RRC = ifft(H_RRC,'symmetric');
h_RRC = ifftshift(h_RRC/sqrt(h_RC(1)));
h_RC = ifftshift(h_RC/h_RC(1)); 
 
dt = 1/Fs;   
tvector = (-(N-1)/2-1:(N-1)/2-1)*dt;

figure;
plot(fvector,fftshift(H_RRC),'r-',fvector,fftshift(H_RRC),'b*');
hold off;
figure("Name","Impulse responses of the RC and RRC filters");
p=zeros(1,4);
p(1)=plot(tvector,h_RC,'r-');hold on;p(2)=plot(tvector,h_RC,'b*');hold on;
p(3)=plot(tvector,h_RRC,'g-');hold on;p(4)=plot(tvector,h_RRC,'g*');
hold off
legend([p(1) p(3)],"RC","RRC");

end

