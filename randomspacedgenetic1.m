close all; clear all
% *********************RANDOMLY SPACED PLANAR ARRAY SYNTHESIS**************
frequency=10;% IN GHz
lambda=300/frequency; %IN mm
k=2*pi/lambda;
array_size=[300,600]; %IN mm
Chrom=1000; %initial population size
iteration=10 %no. of iterations 
res=1; %resolution
mut_rate=15; %MUTATON RATE
%INITIAL POPULATION GENERATION
for chr=1:Chrom
    pts1(:,:,chr)=poissonDiscrandomarray(array_size,.8*lambda,150,0,1.4*lambda);
%     pts1(:,1,chr)=pts1(:,1,chr);
    pts(:,:,chr)=[-pts1(:,1,chr) pts1(:,2,chr);pts1(:,:,chr)];
     pts(:,:,chr)=fliplr( pts(:,:,chr));
end


%*****************GNENETIC ALGORITHM******************************
for iter=1:iteration
   for chr=1:Chrom
u=-1:res/90:1;
    v=-1:res/90:1;
    for i=1:length(u)
        for j=1:length(v) 
            
  AF(i,j,chr)=sum(exp(1j*(2*pi/lambda)*(u(i)*pts(:,1,chr)+v(j)*pts(:,2,chr))));
            end
    end
      AFdB(:,:,chr)=20*log10(abs(AF(:,:,chr))/max(max(abs(AF(:,:,chr)))));   
 

  %CROSSOVER
 %SLL calculation
%elevation cut
lobesel=sort(findpeaks(AFdB(:,1+90/res,chr)),'descend');
SLLel(chr)=lobesel(2);

lobesazi=sort(findpeaks(AFdB(1+90/res,:,chr)),'descend');
SLLazi(chr)=lobesazi(2);

SLerror(chr)=20+0.5*(SLLel(chr)+SLLazi(chr));
     end
%population sorting
  [error idx]=sort(SLerror,'ascend');
  pts(:,1,1:round(0.5*length(idx)))=pts(:,1,idx(1:round(0.5*length(idx))));
  pts(:,2,1:round(0.5*length(idx)))=pts(:,2,idx(1:round(0.5*length(idx))));
     bestofgen(iter)=error(1);
   
   %crossover && mutations
   for kk=2:2:Chrom/2
   pts(:,:,kk-1+round(0.5*length(idx)))=[pts(find(pts(:,1,kk-1)<0 & pts(:,2,kk-1)>0),:,kk-1);pts(find(pts(:,1,kk-1)>0 & pts(:,2,kk-1)<0),:,kk-1);pts(find(pts(:,1,kk-1)>0 & pts(:,2,kk-1)>0),:,kk);pts(find(pts(:,1,kk-1)<0 & pts(:,2,kk-1)<0),:,kk)];
    pts(:,:,kk-1+round(0.5*length(idx)))=[pts(find(pts(:,1,kk)<0 & pts(:,2,kk)>0),:,kk);pts(find(pts(:,1,kk)>0 & pts(:,2,kk)<0),:,kk);pts(find(pts(:,1,kk)>0 & pts(:,2,kk)>0),:,kk-1);pts(find(pts(:,1,kk)<0 & pts(:,2,kk)<0),:,kk-1)]; 
   
   %Mutation
   for jk=1:length(mut_rate)
       pts(jk,:,kk-1+0.5*round(length(idx)))=pts(jk,:,kk-1+0.5*round(length(idx)))+0.1*rand*lambda;
       pts(jk,:,kk+0.5*round(length(idx)))=pts(jk,:,kk+0.5*round(length(idx)))+0.1*rand*lambda;
   
   end 
   end
   iter
end

figure(1)
plot(AFdB(90/res+1,:,idx(1)))

figure(2)
plot(AFdB(:,90/res+1,idx(1)))

figure(3)
 plot(pts(:,1,idx(1)),pts(:,2,idx(1)),'o')
 
 figure(4)
 plot(bestofgen)
    