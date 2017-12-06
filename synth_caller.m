function synthesized_sig = synth_caller(speech_file,pitch_file,max_v_freq_file,frame_map_file,win_size,FFT_size,lpc_poles,ar_poles)
    [x fs] = audioread(speech_file);
    signal=x(:,1);
    
    pitch = load(pitch_file);
    pitch = pitch(:,2);
    load(max_v_freq_file);
    max_v_freq=MAX_VOICED_FREQ;
    
%     maxv=medfilt1(max_v_freq,100);  % sMOOTHING  THE MAX_V_FREQ
    %max_v_freq(max_v_freq>130*20)=130*20;
    
    frame_map=csvread(frame_map_file);
    frame_map=frame_map(:,1);
    shift=win_size/2;
    dest_frames=frame_map;
    [amplitudes,phases,corrected_phi,harm_number]=param_gen(signal,fs,pitch,max_v_freq,frame_map,win_size,FFT_size,lpc_poles);
    synthesized_sig=ola(signal,fs,dest_frames,shift,win_size,pitch,max_v_freq,amplitudes,phases,corrected_phi,harm_number,FFT_size,ar_poles);
%     
%     plot(synthesized_sig);
%     hold on;
%     plot(signal);
end
            
            
            
    