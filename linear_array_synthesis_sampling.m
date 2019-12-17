close all
clear all
% ****************************LINEAR ARRAY SYTHESIS USING SAMPLING********
% Note: Real coefficients can be obtained for symmetric patterns only(Line 42)
N1=21; %NUMBER OF ELEMENTS
frequency=10;% IN GHz
lambda=300/frequency; % IN mm
k=2*pi/lambda;
d=lambda/2; %INTER-ELEMENT SPACING
PatterndB=load('idealpattern.txt');%LOAD INPUT PATTERN(IN DB)
% PatterndB(901:1801)=flipud(PatterndB(1:901));

% ************************************************************************
% yyy=-90:90;       UNCOMMENT IF INPUT PATTERN REQUIRES INTERPOLATION
% yy=PatterndB1(:,2)'-max(PatterndB1(:,2)); 
% xx=PatterndB1(:,1)';
% PatterndB=spline(xx,yy,yyy);
Pattern=10.^(PatterndB/20);
theta=0:180/(length(Pattern)-1):180;
sample=0;
cos_theta_m=0;
g=cosd(theta);
for i=1:2*N1-1 %2*N1-1 is the total number of samples required
   [w I]=min(abs(g-(i+0.5-N1)/(N1-1)));
   cos_theta_m=[cos_theta_m g(I)];
   sample=[sample Pattern(I)];
end
sample=sample(2:length(sample));
cos_theta_m=cos_theta_m(2:length(cos_theta_m));
N=floor(N1/2);
if rem(N1,2)==0
zn=-N:N-1;
else
zn=-N:N;
end

zn=k*zn*lambda/2;
power=zn'*cos_theta_m;
exc=(1/N1)*exp(-1j*power)*sample';
excm=max(abs(exc));
EXCITATION_MAGNITUDE_OUTPUT=abs(exc)/excm
EXCITATION_PHASE_OUTPUT=(180/pi)*acos(real(exc/excm))
REAL_COEFFICIENTS=real(exc); %Only useful if input pattern is even 
z=zeros(N1,length(theta));
     for i1=1:N1
         z(i1,:)=exp(1j*(((N1+1)/2-i1)*(k*d*cosd(theta))));  
     end
  AF=exc'*z;
 AF=abs(AF/max(abs(AF)));
  AF=20*log10(AF);
  a=find(AF<-50);
  AF(a)=-50;
  
  
 figure('name','Array factor plot')
 plot(theta,AF)
hold on
plot(theta,PatterndB)
hold off
legend('synthesized pattern','input pattern')
title('array factor')




