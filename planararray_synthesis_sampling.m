
close all
clear all
%*******************PLANAR ARRAY SYNTHESIS USING SAMPLING*****************************************
% ***************************************INPUTS************************************************************
M=20; 
M1=2*M; %NUMBER OF ROWS 
N=20; 
N1=2*N; %NUMBER OF COLOUMNS
frequency=10;%IN GHz
lambda=300/frequency; % IN mm
k=2*pi/lambda;
PatterndB1=load('Azi.txt');%LOAD AZIMUTH PATTERN(IN dB),THETA IN FIRST COLUMN, DESIRED PATTERN IN 2ND COLUMN
PatterndBel1=load('eleva.txt'); %LOAD ELEVATION PATTERN(IN dB),THETA IN FIRST COLUMN, DESIRED PATTERN IN 2ND COLUMN
delta=lambda/(2*lambda);  %MINIMUM SAMPLING INTERVAL
samplingpointseperation=delta/10000; %SAMPLING INTERVAL
resolution=1; %RESOLUTION
AFres=1; %RESOLUTION OF AF PLOT
dx=lambda/2;%inter element spacing along x
dy=lambda/2;%inter element spacing along y
%********************************************************************
%ELEVATION PATTERN CUBIC INTERPOLATION FOR SAMPLING
X=-90:resolution:90;
yy=PatterndBel1(:,2)'-max(PatterndBel1(:,2));
xx=PatterndBel1(:,1)';
PatterndBel=spline(xx,yy,X);
Patternel=10.^(PatterndBel/20); 
%AZIMUTH PATTERN CUBIC INTERPOLATION FOR SAMPLING
X=-90:resolution:90;
yy=PatterndB1(:,2)'-max(PatterndB1(:,2));
xx=PatterndB1(:,1)';
PatterndB=spline(xx,yy,X);
Pattern=10.^(PatterndB/20);


%SAMPLING IN AZIMUTH
sample=0;
cos_theta_m=0;
theta=0:180/(length(Pattern)-1):180;
g=cosd(theta);
for i=1:floor(2/samplingpointseperation)
   [w I]=min(abs(g-(-1+(i-1)*samplingpointseperation)));
   cos_theta_m=[cos_theta_m g(I)];
   sample=[sample Pattern(I)];
end
sample=sample(2:length(sample));
cos_theta_m=cos_theta_m(2:length(cos_theta_m));
zn=-N+1:N;
zn=k*zn*lambda/2;
power=zn'*cos_theta_m;
exc=(1/N1)*exp(-1j*power)*sample';
excm=max(abs(exc));
exc=exc/excm;

%SAMPLING IN ELEVATION 
sample1=0;
cos_theta_m1=0;
theta1=0:180/(length(Patternel)-1):180;
g1=cosd(theta1);
g=cosd(theta);
for i=1:floor(2/samplingpointseperation) %2*N1-1 is the total number of samples required
   [w I]=min(abs(g-(-1+(i-1)*samplingpointseperation)));
   cos_theta_m1=[cos_theta_m1 g(I)];
   sample1=[sample1 Patternel(I)];
end
sample1=sample1(2:length(sample1));
cos_theta_m1=cos_theta_m1(2:length(cos_theta_m1));
zn1=-M+1:M;
zn1=k*zn1*lambda/2;
power1=zn1'*cos_theta_m1;
excel=(1/M1)*exp(-1j*power1)*sample1';
excmel=max(abs(excel));
excel=excel/excmel;

%EXCITATION COEFFICIENTS FOR PLANAR ARRAY
excoeff=(excel*exc')
excoeffmax=max(max(real(excoeff)));
excoeff=(real(excoeff))/excoeffmax;
MAGNITUDE_OUTPUT=abs(excoeff)%MAGNITUDE OF EXCITATION
PHASE_OUTPUT=(180/pi)*angle(excoeff) %PHASE OF EXCITATION(IN DEGREES)
%excoeff=excoeff(M/4:3*M/4+1,N/4:3*N/4);
%*************** ARRAY FACTOR CALCULATION USING SYNTHESIZED COEFFICIENTS**************************
 u=-1:AFres/90:1;
     v=-1:AFres/90:1;
    
        for f=1:180/AFres+1
            for g=1:180/AFres+1
                x=exp(1j*(k*dx*u(f)));
                y=exp(1j*(k*dy*v(g)));

              Q1y=1*y.^(0:N1-1);
              Qx1=1*x.^(0:M1-1);
              Q=Qx1'*Q1y;

              AF(f,g)=abs(sum(sum(excoeff(:,:).*Q)));

             end
        end
        AF=20*log10(AF/max(max(AF)));
        indices=find(AF<-70); %truncating at -70dB
        AF(indices)=-70;
        
%  ************************PLOTS*************************************8
        figure('Name','2D ARRAY FACTOR PLOT','NumberTitle','off')
        mesh(u,v,AF)
        xlabel('azimuth')
        ylabel('elevation')
        
        figure('Name','Azimuth Cut','NumberTitle','off')
        plot(PatterndB)
        hold on
        plot(AF(1+90/AFres,:))
        hold off
        legend('input pattern','synthesized pattern')
        
        figure('Name','Elevation Cut','NumberTitle','off')
        plot(PatterndBel)
        hold on
        plot(AF(:,1+90/AFres))
        hold off
         legend('input pattern','synthesized pattern')
        