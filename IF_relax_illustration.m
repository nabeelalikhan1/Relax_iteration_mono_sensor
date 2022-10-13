clear all;
close all;
N_sensors=128/32;
N_sensors=4;
n=0:127;
%addpath('D:\D\win64_bin\win64_bin');
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('E:\Published Papers\DOA estimation of intersecting components 2018\Matlab code');
%addpath('E:\Published Papers\DOA ESTIMATION VITERBI\Multi-sensor IF estimation code');

%crossing componentsi8

s1=1.*exp(2*pi*1i*(0.05*n+0.3*n.^3/(128*128*3)));
%s2=1*exp(2*pi*1i*(0.32*n-0*0.3*n.^3/(128*128*3)));
s2=1.*exp(2*pi*1i*(0.09*n+1*0.3*n.^3/(128*128*3)));

s3=1.*exp(2*pi*1i*(0.4*n-0.35*n.^3/(128*128*3)));
%s2=1.*exp(2*pi*1i*(0.35*n+0.3*n.^2/(128*8)));
%s5=1.*exp(2*pi*1i*(0.46*n-1*0.35*n.^3/(128*128*3)));
SampFreq=128;
FFT_len=128;
perc=0.4;
%s = [(s1.') (s2.') (s3.') (s4.') (s5.')];%  (s5.') (s6.') (s7.') ];
%s=real(s);
IF_O(1,:)=0.05+0.3*3*n.^2/(128*128*3);
%IF_O(2,:)=0.35+1*0.3*2*n.^1/(128*8);
IF_O(2,:)=0.09+1*0.3*3*n.^2/(128*128*3);
%IF_O(4,:)=0.46-1*0.35*3*n.^2/(128*128*3);
IF_O(3,:)=0.4-3*0.35*n.^2/(128*3*128);
%IF_O=IF_O.';

s=s1+s2+s3;

n_sources=7-2;
s_orig=s;

LL=500;
index=0;
num=3;

X=awgn(s,2,'measured');

win_length=65;
FFT_len=128;
L=64;
delta=2;
%I=HTFD_new1(Sig(1,:),2,8,64);
 %     figure;
  %   imagesc(I)
    
[ IF1,Xout ] = relax_filtering_mono_sensor( X,num,win_length,delta,L,1,FFT_len );
figure
plot(IF_O.','b');
hold on;
plot(IF1.','r:')
title('Relax IF estimation')



[IF2,Xout] = FASTEST_IF(X,win_length, num, delta,L,0,0,1,FFT_len);
figure
plot(IF_O.','b');
hold on;
plot(IF2.','r:')
title('IF estimation')
