close all
clear all
%***FAILURE CORRECTION FOR LINEAR ARRAY USING ITERATIVE FOURIER TRANSFORM**
%%******INPUTS*********************************************************
M=40; %                                             NUMBER OF ELEMENTS
frequency=10;%                                      FREQUENCY(IN GHz)
lambda=300/frequency; %                               (IN mm)
k=2*pi/lambda;
d=lambda/2;%                                        INTEL ELEMENT SPACING
specSLL=-17;%                                      SPECIFIED SIDELOBE LEVEL
specnulldepth=-100;%                                SPECIFIED NULL DEPTH
minimum_excitation=0.1; %                          MINIMUM EXCITATION AMPLITUDE
failed_elements=load('failed_elements.txt');  %INPUT FAILED ELEMENT INDEX
exc=ones(1,M);
nopoints=8000;%                               n-POINT FOURIER TRANFORM(RESOLUTION)
last=10000;%                                  Number of iterations


% ********************FOURIER TRANSFORM ITERATIONS**********************************
for iter=1:last
    
AF=ifftshift(ifft(exc,nopoints));
AFabs=abs(AF);
AFmax=max(max(AFabs));
AFdB=20*log10(abs(AF)/AFmax);
if iter==1
    AF11=AF;
    AF1=AFdB; %       IDEAL PATTERN
end
if iter==2  %         FAILED PATTERN
    AF22=AF;
    AF2=AFdB;
end
if iter==last-1  %         RESUMED PATTERN
    AF33=AF;
    AF3=AFdB;
end
[nulls indmin]=findpeaks(-AFabs);  %NULLS
[peaks indpeak]=findpeaks(AFabs); %PEAKS
[peaklevel indP]=sort(peaks,'descend');
indpeak=indpeak(indP);

SLL=peaklevel(2:length(peaklevel));
indSLL=indpeak(2:length(indpeak));
specSLL=specSLL;
adapt=find(SLL>max(max(AFabs))*10.^(specSLL/20)); 
AF(indSLL(adapt))=max(max(AFabs))*10.^specSLL/20;% ADAPTING FAILED PATTERN SLL
[depth null_lowering]=findpeaks(-AFabs);
 AF(null_lowering)=max(max(AFabs))*10^(specnulldepth/20);   %    ADAPTING FAILED PATTERN NULLS

exc=fft(ifftshift(AF));    %                     INVERSE FOURIER TRANSFORM
exc=exc(1:M); %                                    truncating excitation coefficients
exc=exc/(max(abs(exc)));
if iter==last-1
    exc_resumed=exc;
end
excM=abs(exc)/max(abs(exc));
excP=(180/pi)*angle(exc);
iter
b=find(excM<minimum_excitation);                               %ADAPTING TO MINIMUM EXCITATION AMPLITUDE
exc(b)=minimum_excitation*exc(b)./abs(exc(b));

exc(failed_elements)=0;                    % SETTING FAILED ELEMENTS TO ZERO;
end


fidRx = fopen('Mag&phase.txt','w');
fprintf(fidRx,'%f\t%f\n',excM,excP);
a=find(AF1<-50);
AF1(a)=-50;
figure
plot(AFdB)
hold on
plot(AF2)
hold off
legend('AF recovered','AF after failure')

figure
plot(AFdB)
hold on
plot(AF1)
legend('AF recovered','AF ideal')
hold off

figure
plot(excM)
