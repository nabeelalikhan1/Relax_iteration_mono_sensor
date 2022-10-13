clear all;
close all;
n=0:127;
s1=1.*exp(2*pi*1i*(0.00*n+0.3*n.^3/(128*128*3)));
%s2=1*exp(2*pi*1i*(0.32*n-0*0.3*n.^3/(128*128*3)));
s2=1.*exp(2*pi*1i*(0.1*n+1*0.3*n.^3/(128*128*3)));

s3=1.*exp(2*pi*1i*(0.4*n-0.35*n.^3/(128*128*3)));
%s2=1.*exp(2*pi*1i*(0.35*n+0.3*n.^2/(128*8)));
%s5=1.*exp(2*pi*1i*(0.46*n-1*0.35*n.^3/(128*128*3)));
%s = [(s1.') (s2.') (s3.') (s4.') (s5.')];%  (s5.') (s6.') (s7.') ];
%s=real(s);
IF_O(1,:)=0.0+0.3*3*n.^2/(128*128*3);
%IF_O(2,:)=0.35+1*0.3*2*n.^1/(128*8);
IF_O(2,:)=0.1+1*0.3*3*n.^2/(128*128*3);
%IF_O(4,:)=0.46-1*0.35*3*n.^2/(128*128*3);
IF_O(3,:)=0.4-3*0.35*n.^2/(128*3*128);
SampFreq=128;
FFT_len=128;
perc=0.4;
win_length=65;
L=32;
%L=5;
s=s1+s2+s3;
LL=250;
index=0;
delta=1;
num=3;
addpath('D:\tfsa_5-5\windows\win64_bin');
for SNR=0:2:10
    %for SNR=-2:2:0
    for ii=1:LL
        
        
        %A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  % mixing matrix A
                                   % mixed source
        % generate noise
        
        X=awgn(s,SNR,'measured');
        
        
        %tic
        for k=0:2
            
            switch k
                case 0
                    [ IFF,~ ] = relax_filtering_mono_sensor( X,num,win_length,delta,L,1,FFT_len );
                    
                case 1
                    
                    %[IFF,ss] = Multi_Sensor_FASTEST_IF_new(X,N_sensors,65, n_sources, delta,64,0,0,4,length(X));
                    [ IFF,~ ] = FASTEST_IF(X,win_length, num, delta,L,0,0,1,FFT_len);
                case 2
                     IFF =ADTFD_RANSAC(X,3,15,64,num,500/1,4,32,8);
                     IFF=IFF/(2*length(X));
                    
            end
            
            
            
            %toc
            
            
            msee=0.1*ones(1,num);
            
            for ii22=1:num
                
                t=1:128;
                IF=IFF(ii22,:);%/length(X);
                t=t(5:end-5);
                for i=1:num
                    c(i)=sum(abs(IF(t)-IF_O(i,t)).^2);
                end
                [a1, b1]=min(c);
                if msee(b1)>=a1(1)/length(X)
                    msee(b1)=a1(1)/length(X);
                end
                
            end
            
            switch k
                case 0
                    mseeIF_refined(ii)=mean(msee);
                case 1
                    mseeIF_new(ii)=mean(msee);
                case 2
                    mseeIF_adtfd1(ii)=mean(msee);
                    
            end
            
        end
    end
    index=index+1;
    %mean(mmssee)
    snr_mse_refined(index)=mean(mseeIF_refined)
    snr_mse_new(index)=mean(mseeIF_new)
    snr_adtfd(index)=mean(mseeIF_adtfd1)
    
    
    
end

SNR=0:2:10;
plot(SNR,10*(log10(snr_mse_refined)),'--md','linewidth',3);
hold on;
plot(SNR,10*(log10(snr_mse_new)),'k','linewidth',3);
hold on;
plot(SNR,10*(log10(snr_adtfd)),'b:','linewidth',3);


xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
legend('Refined IF-estimation using 2-stage method','IF estimation using FAST-IF','ADTFD-RANSAC');
%legend('The Proposed Method','Time-frequency Music','DOA based on IF estimation using ridge tracking');

