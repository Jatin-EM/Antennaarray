close all
clear all

%%INPUTS
N=20; %No. of rows
M=20;%No. of columns
frequency=10;%(GHz)
lambda=300/frequency; %(mm)
k=2*pi/lambda;
dx=lambda/2;%inter element spacing along x
dy=lambda/2;%inter element spacing along y
m=1800; %N-point fft
n=m;
specSLL=-16;
specnulldepth=-160;
exc=ones(M,N);
failed_ele_x=randi([1,M],1,80);
failed_ele_y=randi([1,M],1,80);

for iter=1:100
 
AF=ifftshift(ifft2(exc,m,n));
AFabs=abs(AF);
AFmax=max(max(abs(AF)));
 AFdB=20*log10(abs(AF)/max(max(abs(AF))));
 if iter==1
     AF1=AFdB;
     a1=find(AF1<-70);
     AF1(a1)=-70;
 end
 if iter==2
     AF2=AFdB;
     a1=find(AF2<-70);
     AF2(a1)=-70;
 end
 [pksy idy]=findpeaks(AFabs(round(m/2)+1,:));
 [pksx idx]=findpeaks(AFabs(:,1+round(m/2)));
 [nullx idy1]=findpeaks(-AFabs(1+round(m/2),:));
 [nully idx1]=findpeaks(-AFabs(:,1+round(m/2)));
 [p1 mbx]=sort(pksx,'descend');
 [p2 mby]=sort(pksy,'descend');
 
%  idy=[91*ones(length(idy),1) idy];
if AFdB(1+round(m/2)+2,idy)>specSLL
 AF(1+round(m/2)+2,idy)=max(max(AFabs))*10^(specSLL/20);
end
if  AF(1+1+idx,round(m/2)+2)>specSLL
 AF(1+1+idx,round(m/2)+2)=max(max(AFabs))*10^(specSLL/20);
end
  AF(1+round(m/2),idy1)=max(max(AFabs))*10^(specnulldepth/20);
 AF(idx1,1+round(m/2))=max(max(AFabs))*10^(specnulldepth/20);
 AF(idx(mbx(1)),idy(mby(1)))=AFmax;
 exc=fft2(ifftshift(AF));
%   exc=exc(round(m/2)-25:25+round(m/2),round(m/2)-25:round(m/2)+);
exc=exc(1:M+1,1:N+1);
 excmax=max(max(abs(exc)));
 exc=exc/excmax;
 excM=abs(exc);
 excP=(180/pi)*angle(exc);
 b=find (excM<0.1);
 exc(b)=0.1*exc(b)./excM(b);
 for jj=1:80
 exc(failed_ele_x(jj),failed_ele_y(jj))=0;
 end
 iter
 
end

for jj=1:80
 a(jj)=excM(failed_ele_x(jj),failed_ele_y(jj));
 end
figure
plot(a)

a1=find(AFdB<-70);
     AFdB(a1)=-70;

figure
plot(AFdB(round(m/2)+2,:))
hold on
plot(AF2(round(m/2)+1,:))
hold off
legend('AZIresumed pattern','failed pattern')

figure
plot(AFdB(:,round(m/2)+2))
hold on
plot(AF1(:,round(m/2)+1))
hold off
legend('ELEresumed pattern','ideal pattern')

figure
plot(AFdB(:,round(m/2)+2))
hold on
plot(AF2(:,round(m/2)+1))
hold off
legend('ELEresumed pattern','failed pattern')

figure
plot(AFdB(round(m/2)+2,:))
hold on
plot(AF1(round(m/2)+1,:))
hold off
legend('AZIresumed pattern','ideal pattern')