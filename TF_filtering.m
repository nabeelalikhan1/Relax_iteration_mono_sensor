
function [x,Sig] = TF_filtering(Sig,IF,LL)

Phase=2*pi*filter(1,[1 -1],IF);
    s_dechirp=exp(-1i*Phase);
    
        
                [x,~,~] = ICCD_sparse(Sig,1,IF,3,1,1:length(Sig));

        %wo=asind((ww-1)/64);
        
    %    Sig(jj,:)=s2;%-extr_Sig(iii);
 %       Xout(jj,:)=s1;
        % end
    
    %A = conj(Vector)/norm(Vector);
%A=conj(A);
   Sig=Sig-x;%*abs(sum(AA*X)).^2;
   %I=HTFD_new1(s1,3,8,64);
   %figure; imagesc(I)
  % Sig=Sig-pinv(A)*x;%*abs(sum(AA*X)).^2;
   
  % x=pinv(A.')*Xout1;
   
   %Sig=Sig-A.'*x;
   
   %x = A*Xout1;%/N_S;
  % Sig=Sig-pinv(A)*x;%*abs(sum(AA*X)).^2;
   
   
end



