function [ IF,Xout ] = relax_filtering_mono_sensor( X,n_sources,win_length,delta,L,step,FFT_len )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Sig=X;

for i=1:n_sources
[IF(i,:),Xout] = FASTEST_IF(Sig,win_length, 1, delta,L,0,0,step,FFT_len);
Sig=Sig-Xout;
end
%figure; plot(IF.')
%IF1=IF;
%Veccc1=Veccc;
    for i2=1:3
        % IF=IF1;
         %   Veccc=Veccc1;
           
        for i=1:n_sources
            Sig=X;
            for j=1:n_sources
                
                if i~=j   % REmove all components except i-th
                    [x,~] = TF_filtering(X,IF(j,:),2);
                    Sig=Sig-x;
                    %[Sig] = TF_SF_filtering_new(Sig,Xout(i,:),Veccc(j,:));

                end
                
                
            end
            % Restimate i-th component
              [IF(i,:),Xout(i,:)] = FASTEST_IF(Sig,win_length, 1, delta,L,0,0,step,FFT_len);
%I=HTFD_new1(Xout(i,:),3,8,64);
%figure;imagesc(I)
        end
    end


end

