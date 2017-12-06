function synthesized_sig = ola(signal,fs,dest_frames,shift,win_size,pitch,max_v_freq,amplitudes,phases,corrected_phi,harm_number,FFT_size,ar_poles)
     
     synthesized_sig = zeros(1,dest_frames(end)*shift+win_size);
     BPFILT=designfilt('bandpassfir','FilterOrder', 101, 'CutoffFrequency1', 4000, 'CutoffFrequency2', 7400,'SampleRate', fs);

    for frame_index=1:length(dest_frames)-1
        
        source_frame=dest_frames(frame_index);
        START=(source_frame-1)*shift+1; %Starting point of input frame
        BEGIN=(frame_index-1)*shift+1;  % Starting point of output frame
        
        sig=signal(START:START + win_size-1);
        sig_var=var(sig);
        
        if(pitch(frame_index)==0)
            harm_sig=0;
        else
            %N_harm=floor(max_v_freq(source_frame)/pitch(frame_index));

            a_k=amplitudes{frame_index};
            phi_k=phases{frame_index};
            phi_corr=corrected_phi{frame_index};
            N_harm=harm_number{frame_index};
            harm_sig=harmonic_gen(a_k,phi_k,phi_corr,fs,pitch(frame_index),win_size,FFT_size,N_harm,sig_var);
                  
        end
        [noise_sig,noise_var] = noise_gen(sig,fs,win_size,ar_poles,max_v_freq(source_frame),BPFILT,harm_sig);

        if(pitch(frame_index)~=0)
            if(var(sig)>noise_var)
                
                var_new=var(sig)-noise_var;
                %var_new=var(sig);
            else
                var_new=var(sig);
            end
             harm_sig=harm_sig.*sqrt(var_new);
        end
%         synth=(noise_sig+harm_sig).*hanning(win_size).';
%         synth=(noise_sig+harm_sig).*triang(win_size).';
        synth=(noise_sig+harm_sig);
        
        
        if(frame_index==1)
            synthesized_sig(BEGIN:BEGIN+win_size-1)=synth;
            old_synth=[synth(shift+1:end) zeros(1,shift)];
        else
            synthesized_sig(BEGIN:BEGIN+win_size-1)=old_synth+synth;
            old_synth=[synth(shift+1:end) zeros(1,shift)];
        end
        
    end
%     plot(synthesized_sig);
%     hold on;
%     plot(signal);
end
            
            
            
    