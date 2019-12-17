clear all
close all
%***SYNTHESIS OF LINEAR RANDOMLY SPACED ARRAY USING STATISTICAL DENSITY TAPER*********
%*****************INPUTS***************************
N=22;% NUMBER OF ELEMENTS(SHOULD BE EVEN)
f=10;% FREQUENCY
resolution=0.01;
lambda=300/f;
arraylength=16; %REQUIRED LENGTH OF ARRAY (IN TERMS OF WAVELENGTH)
reqSLL=30;% REQUIRED SIDE LOBE LEVEL
 y=taylorTappEven(N/2,reqSLL,3);
 min_spacing=0.5 % MINIMUM INTER-ELEMENT SPACING IN TERMS OF WAVELENGTH
 max_spacing=1.1 %MAXIMUM INTER-ELEMENT SPACING IN TERMS OF WAVELENGTH
 
 
%  *************************************************************************
% y=abs(exc);
x=1:length(y);
x1=1:resolution:length(y);
 y1=spline(x,y,x1);
 A=trapz(y1);
 a1=A/N;
 for i=1:length(y1)
yy(i)=(trapz(y1(1:i))/a1)-floor(trapz(y1(1:i))/a1);
 end
 [y2 xx]=findpeaks(-yy);

xx=[1 xx(1:end) length(y1)];
for j=1:length(xx)-1
x_pos(j)=0.5*(xx(j)+xx(j+1));
end

figure
plot(x1,y1)
hold on
plot(1+xx(1:end)*resolution,y1(xx(1:end)),'o')
hold on
plot(1+x_pos*resolution,ones(1,length(x_pos)),'*')
hold off
x_pos=arraylength*(x_pos*resolution/N); %normalizing with wavelength
delx=(x_pos)-[0 x_pos(1:end-1)];
for ij=1:length(delx)
   if delx(ij)<min_spacing
       delx(ij)=min_spacing;
   end
   if delx(ij)> max_spacing
           delx(ij)=max_spacing;
   end
end
delx
%%ARRAY FACTOR
theta=0:resolution:180;
z=zeros(N,length(theta));
     for i1=1:N
         z(i1,:)=exp(1j*(2*pi*sum(delx(1:i1))*cosd(theta)));  
     end
     phi=-90:resolution:90;
     elempat=cosd(phi);
exc=ones(1,N);
AF=exc*z;
 AF=abs(AF);
  AF=elempat.*AF;
  AA=max(max(abs(AF)));
  AF=20*log10(AF/AA);
  
 
figure
plot(AF)

finalarraylength=sum(delx)