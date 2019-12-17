close all; clear all
array_size=[400,400]; %IN mm

frequency_high=10;% IN GHz
lambda_high=300/frequency_high; %IN mm
k_high=2*pi/lambda_high;
frequency_low=10;% IN GHz
lambda_low=300/frequency_low; %IN mm
k_low=2*pi/lambda_low;

Chrom=2; %initial population size
iteration=1 %no. of iterations 
res=0.1; %resolution
mut_rate=10; %MUTATON RATE

% intital population generation
for chr=1:Chrom
    pts_low1(:,:,chr)=poissonDiscrandomarray(array_size,0.5*lambda_low,50,0,lambda_low);
    pts_low1(:,1,chr)=0.25*lambda_low+pts_low1(:,1,chr);
    pts_low(:,:,chr)=[-pts_low1(:,1,chr) pts_low1(:,2,chr);pts_low1(:,:,chr)];
end

for chr=1:Chrom
    pts_high1(:,:,chr)=poissonDischighfreq(array_size,0.5*lambda_high,50,0,pts_low);
    pts_high1(:,1,chr)=0.25*lambda_high+pts_high1(:,1,chr);
    pts_high(:,:,chr)=[-pts_high1(:,1,chr) pts_high1(:,2,chr);pts_high1(:,:,chr)];
end

figure(3)
plot(pts_low(:,1),pts_low(:,2),'*')
hold on
plot(pts_high(:,1),pts_high(:,2),'o')
hold off