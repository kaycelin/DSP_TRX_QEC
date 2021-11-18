%% 2020-01-22, LO amplitude normalize
%% 2020-02-18, Issue1, lo normalize will impact the IQ Magnitude Imbalance
%% 2020-02-18, Issue1 fixed, Add RMS_lo
%% 2020-02-25, Export waveform_I and waveform_Q to QEC function
%% 2020-03-01, g3, Apply FIR HPF for Upconversion
%% 2020-03-01, g3, Apply FIR LPF for Downconversion
%% 2020-03-28, DSP_filter_g, Modify filtering form timedomain to FreqDomain
%% 2020-04-13, switch waveform to COLUMN: DIM_FFT = 1
%% 2020-09-16, Seperate Downconversion for I and Q, waveform input is [I_Column; Q_Column]
%% 2021-03-11, Based on 2020-09-16, setup waveform's MAX_Size_Column=2
%% 2021-03-13, Remove 2021-03-11

function [waveform_IQ, waveform_I, waveform_Q] = Mixer_Up_Down_Convert_g3(waveformInput, lo_I, lo_Q, fs, disp, flag_Up_Down, flag_FIR, fnum)
%% Double Balance Mixer

% COLUMN
DIM_FFT = 1;
if size(waveformInput,1)>size(waveformInput,2) % COLUMN
    flag_wf_original = 'COLUMN';
    lo_I=lo_I(:); % switch to COLUMN
    lo_Q=lo_Q(:); % switch to COLUMN
else
    flag_wf_original = 'ROW';
    waveformInput=waveformInput.'; % switch to COLUMN
    lo_I=lo_I(:); % switch to COLUMN
    lo_Q=lo_Q(:); % switch to COLUMN
end

%% 2021-03-13, Remove 2021-03-11
% %% 2021-03-11, Based on 2020-09-16, setup waveform's MAX_Size_Column=2
% if size(waveformInput,DIM_FFT+1)>2
%     error('waveform MAX_Size_Column=2, [I_Column, Q_Column]!')
% end
    
Nsamps = length(waveformInput);

if ~exist('disp','var') || isempty(disp)
    disp = [];
else
    disp = [disp];
end

if ~exist('fnum','var') || isempty(fnum)
    fnum = 0;
end

if ~exist('lo_I','var') || isempty(lo_I)
    lo_I = 1;
end

if ~exist('lo_Q','var') || isempty(lo_Q)
    lo_Q = 1;
end

%% 2020-02-18, Issue1, lo normalize will impact the IQ Magnitude Imbalance
flag_LO_MAG_Normalize = 'YES'; % default to avoid the Mixer conversion gain impact the pwr saturation

switch flag_LO_MAG_Normalize
    case {'YES'}
        %% 2020-02-18, Issue1 fixed, Add RMS_lo
        RMS_lo = mean([sqrt(mean(abs(lo_I).^2)),sqrt(mean(abs(lo_Q).^2))]);
        lo_I = lo_I/RMS_lo;
        lo_Q = lo_Q/RMS_lo;
        disp_LO_rms = [', /RMS(LO)'];
    case {'NO'}
        disp_LO_rms = [', /1'];
end

%% Upconversion
switch flag_Up_Down
    case {'Up'}
        
        lo_I_Up = lo_I;
        lo_Q_Up = lo_Q;
        
        %         if size(waveform,mod(DIM_FFT,2)+1)==1
        if ~isreal(waveformInput)
            waveform_I_Up = real(waveformInput).*lo_I_Up;
            waveform_Q_Up = imag(waveformInput).*lo_Q_Up;
        elseif size(waveformInput,mod(DIM_FFT,2)+1)==2 % [I_Column; Q_Column]
            waveform_I_Up = waveformInput(:,1).*lo_I_Up;
            waveform_Q_Up = waveformInput(:,2).*lo_Q_Up;
        elseif size(waveformInput,mod(DIM_FFT,2)+1)==1
            display('Signle Waveform is NOT IQ!!!')
            waveform_I_Up = waveformInput.*lo_I_Up;
            waveform_Q_Up = waveformInput.*lo_Q_Up;
%             warning('Format of waveform check!')
        end
        
        % pwr check
        bw_inband = [];
        fnum_debug=[];
        [Ipwr_dB_WF_I_Up, ind_inband, Opwr_dB] = Pwr_Inband_g(fft(waveform_I_Up,[],DIM_FFT), fs, bw_inband, 0e6, 'full', fnum_debug);
        [Ipwr_dB_WF_Q_Up, ind_inband, Opwr_dB] = Pwr_Inband_g(fft(waveform_Q_Up,[],DIM_FFT), fs, bw_inband, 0e6, 'full', fnum_debug);
        
        % IQ combination for Upconversion
        % default 1: Sign_IQ_Lo_Phs_Up: I is 90deg behind of Q;
        % default 2: Sign_IQcomb_Up: I+Q
        Sign_IQ_Lo_Phs_Up = -1; % 1: I in 90deg in front of Q; -1: I is 90deg behind of Q
        Sign_IQcomb_Up = +Sign_IQ_Lo_Phs_Up; % -1: I-Q; +1: I+Q
        
        waveform_IQ_Up = waveform_I_Up+(Sign_IQcomb_Up)*waveform_Q_Up;
        
        flag_debug_plot = 0;
        if flag_debug_plot==1
            fnum_debug=str2double(['0312',num2str(fnum)])
            PLOT_FFT_dB_g(waveform_I_Up, fs, Nsamps, ['waveform I'], 'df', 'full', 'pwr', fnum_debug);
            PLOT_FFT_dB_g(waveform_Q_Up, fs, Nsamps, ['waveform Q'], 'df', 'full', 'pwr', fnum_debug);
            PLOT_FFT_dB_g(waveform_IQ_Up, fs, Nsamps, ['waveform IQ'], 'df', 'full', 'pwr', fnum_debug);
            title('Upconversion I and Q')
        end
        
        %% 2020-03-01, g3, Apply FIR HPF for Upconversion
        if ~exist('flag_FIR','var')||isempty(flag_FIR)||any(flag_FIR==0)
            %% 2020-02-25, Export waveform_I and waveform_Q to QEC function
            waveform_IQ = waveform_IQ_Up;
            waveform_I = waveform_I_Up;
            waveform_Q = waveform_Q_Up;
            
            %             % plot
            disp_IQ_Up = [disp, disp_LO_rms];
            %             figure(str2double(['2020070101']))
            %             PLOT_FFT_dB_g(waveform_IQ, fs, Nsamps, disp_IQ_Up, 1, 'full'); hold on
            %             title('Upconversion')
        else
            bw_inband = flag_FIR;
            FIR_f_cuttoff_tolerance = 0.5e6; % < fs/2
            FIR_Norder=10001;
            if size(bw_inband,2)==1
                FIR_Ftype='HPF';
                f_cutoff_L = 0;
                f_cutoff_H = bw_inband(1)-FIR_f_cuttoff_tolerance;
            elseif size(bw_inband,2)==2
                FIR_Ftype='BPF';
                f_cutoff_L = bw_inband(1)-FIR_f_cuttoff_tolerance;
                f_cutoff_H = bw_inband(2)+FIR_f_cuttoff_tolerance;
                %             elseif strcmp(FIR_Ftype,'BSF')
                %                 bw_inband_org = bw_inband;
                %                 bw_image = bw_inband - 2*fLOU;
                %                 bw_inband = bw_image;
                %                 f_cutoff_L = bw_inband(1)-FIR_f_cuttoff_tolerance;
                %                 f_cutoff_H = bw_inband(2)+FIR_f_cuttoff_tolerance;
                %                 bw_inband = bw_inband_org;
            end
            FIR_Wtype='HAN';
            
            bw_fir_ratio_L = 2*f_cutoff_L/fs;
            FIR_w_cutoff_L = bw_fir_ratio_L*pi;
            bw_fir_ratio_H = 2*f_cutoff_H/fs;
            FIR_w_cutoff_H = bw_fir_ratio_H*pi;
            table_input_FIR = table(FIR_Norder, FIR_Ftype, FIR_Wtype, f_cutoff_L, f_cutoff_H);
            
            % FIR generator
            b = fir_windowdesign_g(FIR_Norder, FIR_Ftype, FIR_w_cutoff_L, FIR_w_cutoff_H, FIR_Wtype);
            
            % FIR plot
            if flag_debug_plot~=0
                yyaxis right
                fnum_debug=str2double(['0312',num2str(fnum)])
                PLOT_FFT_dB_g(b*Nsamps, fs, Nsamps, [FIR_Ftype], 'df', 'full', 'pwr', fnum_debug);
                yyaxis left
            end
            
            % apply FIR filtering
            %             waveform_I_Up_HPF = conv(waveform_I_Up,(b/2^0),'same'); % apply filter in time domain
            %             waveform_Q_Up_HPF = conv(waveform_Q_Up,(b/2^0),'same'); % apply filter in time domain
            waveform_I_Up_HPF = DSP_filter_g(b, waveform_I_Up, 'FD');
            waveform_Q_Up_HPF = DSP_filter_g(b, waveform_Q_Up, 'FD');
            
            if flag_debug_plot~=0
                fnum_debug=str2double(['0301',num2str(fix(fs/1e6))])
                PLOT_FFT_dB_g(waveform_I_Up_HPF, fs, Nsamps, ['waveform I wi ',FIR_Ftype], 'df', 'full', 'pwr', fnum_debug);
                PLOT_FFT_dB_g(waveform_Q_Up_HPF, fs, Nsamps, ['waveform Q wi ',FIR_Ftype], 'df', 'full', 'pwr', fnum_debug);
            end
            
            % IQ combination for Upconversion
            % default 1: Sign_IQ_Lo_Phs_Up: I is 90deg behind of Q;
            % default 2: Sign_IQcomb_Up: I+Q
            Sign_IQ_Lo_Phs_Up = -1; % 1: I in 90deg in front of Q; -1: I is 90deg behind of Q
            Sign_IQcomb_Up = +Sign_IQ_Lo_Phs_Up; % -1: I-Q; +1: I+Q
            
            waveform_IQ_Up_HPF = waveform_I_Up_HPF+(Sign_IQcomb_Up)*waveform_Q_Up_HPF;
            if flag_debug_plot~=0
                fnum_debug = str2double(['0301',num2str(fix(fs/1e6))])
                PLOT_FFT_dB_g(waveform_IQ_Up_HPF, fs, Nsamps, ['txwaveform IQ combination wi ',FIR_Ftype], 'df', 'full', 'pwr', fnum_debug);
            end
            
            % export waveform
            waveform_IQ = waveform_IQ_Up_HPF;
            waveform_I = waveform_I_Up_HPF;
            waveform_Q = waveform_Q_Up_HPF;
            disp_IQ_Up = [disp, disp_LO_rms, ' wi ', FIR_Ftype];
            
            %             % plot
            %             if fnum ~=0
            %                 disp_IQ_Down = [disp, disp_LO_rms, ' wi ', FIR_Ftype];
            %                 figure(str2double(['03120',num2str(fnum)]))
            %                 PLOT_FFT_dB_g(waveform_IQ, fs, Nsamps, disp_IQ_Down, 1, 'full'); hold on
            %                 title('Upconversion')
            %             end
        end
        
        % plot
        if fnum~=0
            %             disp_IQ_Down = [disp, disp_LO_rms, ' wi ', FIR_Ftype];
            %             figure(str2double(['04140',num2str(fnum)]))
            %             figure(fnum)
            PLOT_FFT_dB_g(waveform_IQ, fs, Nsamps, disp_IQ_Up, 'df', 'full', 'pwr', fnum);
            title('Upconversion')
        end
        %% Dpconversion
    case {'Down'}
        
        lo_I_Down = lo_I;
        lo_Q_Down = lo_Q;
        
        %         waveform_I_Down = waveform.*lo_I_Down;
        %         waveform_Q_Down = waveform.*lo_Q_Down;
        %         if size(waveform,mod(DIM_FFT,2)+1)==1
%         if isreal(waveform)
%             waveform_I_Down = waveform.*lo_I_Down;
%             waveform_Q_Down = waveform.*lo_Q_Down;
        if size(waveformInput,mod(DIM_FFT,2)+1)==2
            %% 202-09-16, Seperate Downconversion for I and Q, waveform input is [I_Column, Q_Column]
            waveform_I_Down = waveformInput(:,1).*lo_I_Down;
            waveform_Q_Down = waveformInput(:,2).*lo_Q_Down;
        else
            waveform_I_Down = waveformInput.*lo_I_Down;
            waveform_Q_Down = waveformInput.*lo_Q_Down;
%             error('Format of waveform check!')
        end
        
        % pwr check
        bw_inband = [];
        fnum_debug = [];
        [Ipwr_dB_WF_I_Down, ind_inband, Opwr_dB] = Pwr_Inband_g(fft(waveform_I_Down,[],DIM_FFT), fs, bw_inband, 0e6, 'full', fnum_debug);
        [Ipwr_dB_WF_Q_Down, ind_inband, Opwr_dB] = Pwr_Inband_g(fft(waveform_Q_Down,[],DIM_FFT), fs, bw_inband, 0e6, 'full', fnum_debug);
        
        
        % IQ combination for Downconversion
        % default 1: Sign_IQ_Lo_Phs_Down: I is 90deg behind of Q;
        % default 2: Sign_IQcomb_Down: I+jQ
        Sign_IQ_Lo_Phs_Down = +1; % 1: I in 90deg in front of Q; -1: I is 90deg behind of Q
        Sign_IQcomb_Down = Sign_IQ_Lo_Phs_Down; % -1: I-jQ; +1: I+jQ
        
        waveform_IQ_Down = waveform_I_Down-(Sign_IQcomb_Down)*1i*waveform_Q_Down;
        flag_debug_plot = 0;
        if flag_debug_plot~=0
            PLOT_FFT_dB_g(waveform_IQ_Down, fs, Nsamps, [], 'df', 'full', 'pwr', 20200712);
            PLOT_FFT_dB_g(waveform_Q_Down, fs, Nsamps, [], 'df', 'full', 'pwr', 20200712);
            PLOT_FFT_dB_g(waveform_I_Down, fs, Nsamps, [], 'df', 'full', 'pwr', 20200712);

        end
        
        %% 2020-03-01, g3, Apply FIR LPF for Downconversion
        if ~exist('flag_FIR','var')||isempty(flag_FIR)||any(flag_FIR==0)
            %% 2020-02-25, Export waveform_I and waveform_Q to QEC function
            waveform_IQ = waveform_IQ_Down;
            waveform_I = waveform_I_Down;
            waveform_Q = waveform_Q_Down;
            disp_IQ_Down = [disp, disp_LO_rms];
            
            %             fnum = 0;
            
        else
            bw_inband = flag_FIR;
            FIR_f_cuttoff_tolerance = 0.5e6; % < fs/2
            FIR_Norder=10001;
            if size(bw_inband,2)==1
                FIR_Ftype='LPF';
                f_cutoff_L = 0;
                f_cutoff_H = bw_inband(1)-FIR_f_cuttoff_tolerance;;
                %             elseif strcmp(FIR_Ftype,'BSF')
                %                 bw_inband_org = bw_inband;
                %                 bw_image = bw_inband - 2*fLOU;
                %                 bw_inband = bw_image;
                %                 f_cutoff_L = bw_inband(1)-FIR_f_cuttoff_tolerance;
                %                 f_cutoff_H = bw_inband(2)+FIR_f_cuttoff_tolerance;
                %                 bw_inband = bw_inband_org;
            elseif size(bw_inband,2)==2
                FIR_Ftype='BPF';
                f_cutoff_L = bw_inband(1)-FIR_f_cuttoff_tolerance;
                f_cutoff_H = bw_inband(2)+FIR_f_cuttoff_tolerance;
                %             elseif strcmp(FIR_Ftype,'LPF')&size(bw_inband,2)==1
                %                 f_cutoff_L = bw_inband(1)+FIR_f_cuttoff_tolerance;
                %                 f_cutoff_H = 0;
            end
            FIR_Wtype='HAN';
            
            bw_fir_ratio_L = 2*f_cutoff_L/fs;
            FIR_w_cutoff_L = bw_fir_ratio_L*pi;
            bw_fir_ratio_H = 2*f_cutoff_H/fs;
            FIR_w_cutoff_H = bw_fir_ratio_H*pi;
            table_input_FIR = table(FIR_Norder, FIR_Ftype, FIR_Wtype, f_cutoff_L, f_cutoff_H);
            
            % FIR generator
            b = fir_windowdesign_g(FIR_Norder, FIR_Ftype, FIR_w_cutoff_L, FIR_w_cutoff_H, FIR_Wtype);
            
            flag_debug_plot = 0;
            % FIR plot
            if flag_debug_plot~=0
                yyaxis right
                fnum_debug=str2double(['0312',num2str(fnum)])
                PLOT_FFT_dB_g(b*Nsamps, fs, Nsamps, [FIR_Ftype], 'df', 'full', 'pwr', fnum_debug);
                PLOT_FFT_dB_g(b*Nsamps, fs, Nsamps, [FIR_Ftype], 'df', 'full', 'pwr', 101);
                yyaxis left
            end
            
            % apply FIR filtering
            %             waveform_I_Down_LPF = conv(waveform_I_Down,(b/2^0),'same'); % apply filter in time domain
            %             waveform_Q_Down_LPF = conv(waveform_Q_Down,(b/2^0),'same'); % apply filter in time domain
            waveform_I_Down_LPF = DSP_filter_g(b, waveform_I_Down, 'FD');
            waveform_Q_Down_LPF = DSP_filter_g(b, waveform_Q_Down, 'FD');
            
            if flag_debug_plot==1
                fnum_debug=str2double(['0301',num2str(fix(fs/1e6))])
                PLOT_FFT_dB_g(waveform_I_Down_LPF, fs, Nsamps, ['waveform I wi ',FIR_Ftype], 'df', 'full', 'pwr', fnum_debug);
                PLOT_FFT_dB_g(waveform_Q_Down_LPF, fs, Nsamps, ['waveform Q wi ',FIR_Ftype], 'df', 'full', 'pwr', fnum_debug);
            end
            
            % IQ combination for Upconversion
            % default 1: Sign_IQ_Lo_Phs_Up: I is 90deg behind of Q;
            % default 2: Sign_IQcomb_Up: I+Q
            Sign_IQ_Lo_Phs_Up = +1; % 1: I in 90deg in front of Q; -1: I is 90deg behind of Q
            Sign_IQcomb_Up = +Sign_IQ_Lo_Phs_Up; % -1: I-Q; +1: I+Q
            
            waveform_IQ_Down_LPF = waveform_I_Down_LPF-(Sign_IQcomb_Down)*1i*waveform_Q_Down_LPF;
            if flag_debug_plot~=0
                fnum_debug=str2double(['0312',num2str(fnum)])
                PLOT_FFT_dB_g(waveform_IQ_Down_LPF, fs, Nsamps, ['rxwaveform IQ combination wi ',FIR_Ftype], 'df', 'full', 'pwr', fnum_debug);
            end
            
            % export waveform
            waveform_IQ = waveform_IQ_Down_LPF;
            waveform_I = waveform_I_Down_LPF;
            waveform_Q = waveform_Q_Down_LPF;
            disp_IQ_Down = [disp, disp_LO_rms, ' wi ',FIR_Ftype];
            
            %             % plot
            %             if fnum ~=0
            %                 figure(str2double(['03120',num2str(fnum)]))
            %                 disp_IQ_Down = [disp, disp_LO_rms, ' wi ',FIR_Ftype];
            %                 PLOT_FFT_dB_g(waveform_IQ, fs, Nsamps, disp_IQ_Down, 1, 'full'); hold on
            %                 title('Downconversion')
            %             end
        end
        % plot
        if fnum ~=0
            if size(waveformInput,DIM_FFT+1)==2 % [I_Column; Q_Column]
                PLOT_FFT_dB_g([waveform_I], fs, Nsamps, ['I-',disp_IQ_Down], 'df', 'full', 'pwr', fnum);
                PLOT_FFT_dB_g([waveform_Q], fs, Nsamps, ['Q-',disp_IQ_Down], 'df', 'full', 'pwr', fnum);
            else
                PLOT_FFT_dB_g(waveform_IQ, fs, Nsamps, disp_IQ_Down, 'df', 'full', 'pwr', fnum);
            end
            title('Downconversion')
        end
end

% export and recover
if strcmp(flag_wf_original,'ROW')
    waveform_IQ = waveform_IQ.';
    waveform_I = waveform_I.';
    waveform_Q = waveform_Q.';
    
end