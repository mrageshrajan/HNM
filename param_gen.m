function [a_k, phases, phi_k, Num_harm] = param_gen(speech,fs,pitch,max_v_freq,framemap,win_size,FFT_size,lpc_poles)
    

%     [x fs] = audioread(speech_file);
%     signal=x(:,1);
%     
%     pitch = load(pitch_file);
%     frame_map=csvread(frame_map_file)
    signal=speech;
    shift=win_size/2;
    dest_frames=framemap;
    a_k={1,length(dest_frames)};Num_harm={1,length(dest_frames)};
    phases={1,length(dest_frames)};   
    first_voiced_frame=1;    
    
    for frame_index=1:length(dest_frames)
        source_frame=dest_frames(frame_index);
        START=(source_frame-1)*shift+1; %Starting point of input frame
        BEGIN=(frame_index-1)*shift+1;  % Starting point of output frame
        
        
        
        if(pitch(frame_index)~=0)                  % If pitch values are not zero consider  those frames
    
            sig=signal(START:START+win_size-1);
            SF=fft(sig,FFT_size); % DFT using FFT
            SF1=abs(SF(1:FFT_size/2+1)); % Magnitude of DFT(Should take 1 extra sample after the mid point)
            SF2=angle(SF(1:FFT_size/2+1));
            N_harm=round(max_v_freq(source_frame)/pitch(frame_index)); % NO. of harmonics to be considered
            
            if(N_harm)<1
                N_harm=1;
            elseif (N_harm*pitch(frame_index)>=(fs/2)-200)
                N_harm=floor(((fs/2)-200)/pitch(frame_index));
            end
            
            
            harmos= round([pitch(frame_index):pitch(frame_index):N_harm*pitch(frame_index)]*FFT_size/fs);
            phases{frame_index}=SF2(harmos);
         %%  Smoothed Spectral estimates from LPC Spectrum
         
            lp=lpc(sig,26); % To find LPC coefficients
            [H W]=freqz(1,lp,[0:2*pi/FFT_size:pi]);   % To find LPC spectrum: Since LPC is an all-pole model, 'B' coefficient is 1
            %H=SF1.';
            amp_esti=abs(H(harmos)).';
            a_k{frame_index}=amp_esti./(FFT_size/2);
            Num_harm{frame_index}=N_harm;
        else 
            
            a_k{frame_index}=0;
            Num_harm{frame_index}=0;
            phases{frame_index}=0;
            
        end
        
    end
       
   %%phase correction values 
        
        phi_k={1,length(dest_frames)};
        first_voiced_frame=1;
        %phi_k{1}=0;
        
        for frame_index=1:length(dest_frames)
            source_frame=dest_frames(frame_index);
            START=(source_frame-1)*shift+1; %Starting point of input frame
            
            
            
            if(pitch(frame_index)~=0) % If current frame is voiced
            
                N_harm = length(a_k{frame_index});    % Number of harmonics for current frame
                
                if(first_voiced_frame==1)   % if the current frame is the first voiced frame after an unvoiced frame
                    first_voiced_frame=0;   % If the first voiced frame is found, reset the flag
                    %phi_old=0;  % No phase correction required for current frame
                    %phi_old=2*pi*(rand(1,1));
                    phi_old=phases{frame_index};
                    
                    phi_k{frame_index}=phi_old;
                    
                else    %If there is atleast one voiced frame preceding current frame
                    
                    if(length(a_k{frame_index})>length(a_k{frame_index-1}))

                        phi=zeros(size(a_k{frame_index}));
                        phi(1:length(a_k{frame_index-1}))=phi_old+[1:Num_harm{frame_index-1}].'*(2*pi*pitch(frame_index-1)*(win_size/(2*fs)));
                        phi_k{frame_index}=phi-(floor(phi./(2*pi))*2*pi); % phase wrapping
                        %phi_k{frame_index}=phi;
                        phi_old=phi_k{frame_index};

                    else
                        
                        phi_new=phi_old+[1:Num_harm{frame_index-1}].'*(2*pi*pitch(frame_index-1)*(win_size/(2*fs)));
                        phi=phi_new(1:length(a_k{frame_index}));
                        phi_k{frame_index}=phi-(floor(phi./(2*pi))*2*pi); % phase wrapping
                        %phi_k{frame_index}=phi;
                        phi_old=phi_k{frame_index};
                    end
                   
                end
                 
                
            else% Unvoiced frame
                phi_k{frame_index}=0;
                phi_old=0;
                first_voiced_frame=1; % when an unvoiced frame is foubnd, set the flag to 1 (for the next 1st voiced frame)
                
            end
        end
% return {a_k,phi_k};        
        
end

            
                
        
        
    
            