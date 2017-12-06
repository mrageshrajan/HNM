function [noise_sig,noise_var] = noise_gen(sig,fs,win_size,ar_poles,max_v_freq,BPFILT,harm_sig)
    
    if(max_v_freq~=0)   % If the frame is voiced
        
%% Normal HNM Method By Stylianno    
%%
        filtered=filter(BPFILT,sig);
        noise_sig1=filter(BPFILT,randn(1,length(sig))); % White Gaussian noise
        [ar vr]=arcov(filtered,ar_poles);   % Short term AR modeling of the residue signal
        noise_sig = filter(1,ar,noise_sig1);
        noise_var=var(filtered);
        noise_sig=noise_sig./sqrt(var(noise_sig)).*sqrt(noise_var);
        
%% Overlapping Harmonic and Noise Spectrum
%%        
%        [ar_long vr_long]=arcov(filtered,18);  % Long term AR of filtered ORIGINAL speech--experimentally changed the value 
%        filtered2=filter(BPFILT,harm_sig); % Filter the harmonic signal
%        residue_sig=filter(ar_long,1,filtered);    % finding residue from  filtered ORIGINAL signal
%        [ar vr]=arcov(residue_sig,ar_poles);   % Short term AR modeling of the residue signal
        
%        noise_sig = filter(1,ar,noise_sig1);    % AR filtering 


%        noise_sig=noise_sig./sqrt(var(noise_sig)).*sqrt(var(residue_sig));
%        noise_var=var(residue_sig);

%% Another method for variance calculation in Overlapping spectrum case        
%%      
%         filtered=filter(BPFILT,sig);
%         SF2=fft(filtered.*hanning(length(sig)),4096);
%         filtered2=filter(BPFILT,harm_sig); % Filter the harmonic signal
%         SF3=fft((filtered2.').*hanning(length(harm_sig)),4096)/4096;
%         absf2=abs(SF2);
%         absf3=abs(SF3);
%         angsf2=angle(SF2);
%         angsf3=angle(SF3);
%         t=[0:length(sig)-1]/fs;
%         phii=repmat(angsf2,1,length(t))+2*pi*repmat([0:4095].',1,length(t)).*repmat(t,4096,1);
%         absnew=abs(absf2-absf3*16);
%         a=sum(repmat(absnew(2:2048),1,length(t)).*cos(phii(2:2048,:)));
%         a=(2*a+absnew(1)+absnew(2049))/4096;
%         noise_var=var(a);
        
%% Yet another method        
%%        
%         aa=ones(size(abssf));
%         [yy zz]=findpeaks(abssf,'MinPeakHeight',0.002);
%         aa(zz)=0;
%         
%         aa(zz(2:end)-1)=0;
%         aa(zz(1:end-1)+1)=0;
%         aa(zz(3:end)-2)=0;
%         aa(zz(1:end-2)+2)=0;
%         
%         new=real(SF2).*aa+imag(SF2);
%         iff=ifft(new(1:2049),4096,'symmetric');
%         newsig=iff(1:length(sig));
%         noise_var=var(newsig-filtered);
        
        
        
% 
%         new=(realsf-realsf2)+imag(SF2);
%         iff=ifft(new(1:2049),4096,'symmetric');
%         [ar vr]=arcov(iff,ar_poles);
%         noise_sig1=filter(BPFILT,randn(1,length(sig))); % White Gaussian noise
%         noise_sig = filter(1,ar,noise_sig1);    % 
        
%         noise_var=var(iff(1:length(sig)));
%         noise_sig=noise_sig./sqrt(var(noise_sig)).*sqrt(noise_var);
        
    else
        
        [ar vr]=arcov(sig,ar_poles);
        w = randn(1,length(sig)); 
        noise_sig = filter(1,ar,w);
        noisepart=noise_sig.*(sqrt(vr)/sqrt(var(w)));
        noise_sig=noisepart./sqrt(var(noisepart)).*sqrt(var(sig));
        noise_var=var(noise_sig);
    end
    
end
        
       