function harm_sig=harmonic_gen(a_k,phi_k,phi_corr,fs,pitch,win_size,FFT_size,N_harm,sig_var)
    t=(0:win_size-1)/fs;
% % %     if(phi_corr==0)
% % %         %phi=(repmat([1:N_harm].',1,length(t)).*(2*pi*pitch).*repmat(t,N_harm,1));
% % %         rand_phase=2*pi*(rand(N_harm,1)); % Adding a random phase to avoid accidental peaking of sinusoids while summing. 
% % %         phi=repmat(rand_phase,1,length(t))+(repmat([1:N_harm].',1,length(t)).*(2*pi*pitch).*repmat(t,N_harm,1));
% % %   %%%%%      phi=repmat(phi_k,1,length(t))+(repmat([1:N_harm].',1,length(t)).*(2*pi*pitch).*repmat(t,N_harm,1));
% % %  
% % %     else
% % %       %%%%%  phi=repmat(phi_corr+phi_k,1,length(t))+(repmat([1:N_harm].',1,length(t)).*(2*pi*pitch).*repmat(t,N_harm,1));
% % %         phi=repmat(phi_corr,1,length(t))+(repmat([1:N_harm].',1,length(t)).*(2*pi*pitch).*repmat(t,N_harm,1));
% % %     end
    if(length(phi_corr)==1)
        phi=phi_corr+(repmat([1:N_harm].',1,length(t)).*(2*pi*pitch).*repmat(t,N_harm,1));
    else
        phi=repmat(phi_corr,1,length(t))+(repmat([1:N_harm].',1,length(t)).*(2*pi*pitch).*repmat(t,N_harm,1));
    end
    
    %phi=phi-(floor(phi./(2*pi))*2*pi); % phase wrapping
    h=sum(repmat(a_k,1,length(t)).*cos(phi),1);
    %h=sum(repmat(ones(size(a_k)),1,length(t)).*cos(phi),1);
    %h=repmat(a_k(10),1,length(t)).*cos(phi(10,:));
    harm_sig=h./(sqrt(var(h))); %Normalization
    %harm_sig=h;
    %return harm_sig
end
    