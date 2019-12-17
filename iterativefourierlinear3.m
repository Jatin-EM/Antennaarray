close all
clear all

%%INPUTS
M=40; %No. of elements
frequency=10;%(GHz)
lambda=300/frequency; %(mm)
k=2*pi/lambda;
d=lambda/2;%inter element spacing along x
specSLL=-15;
% yyy=-90:90;
% yy=PatterndB1(:,2)'-max(PatterndB1(:,2));
% xx=PatterndB1(:,1)';
% PatterndB=spline(xx,yy,yyy);
% inpatdB=load('idealpattern.txt');
% inpat=10.^inpatdB/20;
% inpat=inpat';

exc=ones(1,M);
% exc=exc.*taylorTappEven(20,30,3)
 nopoints=8000;


for iter=1:10000
AF=ifftshift(ifft(exc,nopoints));
AFabs=abs(AF);
AFmax=max(max(AFabs));
AFdB=20*log10(abs(AF)/AFmax);
if iter==1
    AF1=AFdB;
end
if iter==2
    AF2=AFdB;
end
%derAF=sign(diff(AFabs));
% Nulls
[nulls indmin]=findpeaks(-AFabs); 
%peaks
[peaks indpeak]=findpeaks(AFabs); 
[peaklevel indP]=sort(peaks,'descend');
indpeak=indpeak(indP);
SLL=peaklevel(2:length(peaklevel));
indSLL=indpeak(2:length(indpeak));
adapt=find(SLL>max(max(AFabs))*10.^(specSLL/20));
AF(indSLL(adapt))=max(max(AFabs))*10.^specSLL/20;
[depth null_lowering]=findpeaks(-AFabs);
 AF(null_lowering)=max(max(AFabs))*10^-5;
%SYNTHESIS
% error=abs(inpat-AFabs);
% adapt=find(error>2);
% AF(adapt)=inpat(adapt);
%  AF=10.^(AFdB/20);
exc=fft(ifftshift(AF));
exc=exc(1:M); %truncating excitation coefficients
exc=exc/(max(abs(exc)));
% exc=exc((m-M)/2:(m+M)/2,(m-M)/2:(m+M)/2);
excM=abs(exc)/max(abs(exc));
excP=(180/pi)*angle(exc);
iter
%b=find(excM<0.5);
% if sum(sum(abs(exc(failed_ele_x,failed_ele_y))))<0.4 && length(b)<15
% break
% end
%exc(b)=0.5*exc(b)./abs(exc(b));
%  b=find(excM<0.5)
%  exc(b)=0.8;
exc([4,7,13,16,17,34,25,38,29,31])=0;
end

% avgideal=sum(AF1)/length(AF1);
% avgsynt=sum(AFdB)/;
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
