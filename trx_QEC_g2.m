%% 2020-02-28, Tolearnce Error Detection --> %% 2020-03-12, Ignore !, Due to the filter has been applied.
%% 2020-03-01, IQ waveform QEC correct
%% 2020-03-04, IQ waveform QEC correct for DL part ?
%% 2020-03-16, Modify the filtering before DownConversion
%% 2020-03-28, DSP_filter_g, Modify filtering form timedomain to FreqDomain
%% 2020-04-13, switch waveform to ROW: DIM_FFT = 2



function [IMB_MagDB, IMB_PhsDeg, waveform_IQ_QEC] = trx_QEC_g3(waveform, lo_I_IMB_Down, lo_Q_IMB_Down, fs, bw_inband, flag_QEC_TX_RX, flag_debug_plot)
% flag_QEC_TX_RX = 'TXQEC'/'TXLOQEC'/'RXQEC'

Nsamps = length(waveform);
Nbr = min(size(waveform));

% switch waveform to ROW
DIM_FFT = 2;
if size(waveform,1)>size(waveform,2) % COLUMN
    flag_wf_original = 'COLUMN';
    waveform=waveform.'; % switch to ROW
else
    flag_wf_original = 'ROW';
end
% switch lo to ROW
if size(lo_I_IMB_Down,1)>size(lo_I_IMB_Down,2) % COLUMN
    lo_I_IMB_Down=lo_I_IMB_Down.'; % switch to ROW
    lo_Q_IMB_Down=lo_Q_IMB_Down.'; % switch to ROW
end
%     waveform=waveform(:); % switch to COLUMN
% end

if ~exist('flag_debug_plot','var')||isempty(flag_debug_plot)
    flag_debug_plot = 0;
end

if ~exist('flag_QEC_TX_RX','var')||isempty(flag_QEC_TX_RX)||strcmp(flag_QEC_TX_RX,'TXQEC')||strcmp(flag_QEC_TX_RX,'TXLOQEC')
    %     if  ~exist('flag_QEC_TX_RX','var')||isempty(flag_QEC_TX_RX)||strcmp(flag_QEC_TX_RX,'TXQEC')
    %         flag_QEC_TX_RX = 'TXQEC';
    %     elseif strcmp(flag_QEC_TX_RX,'TXLOQEC')
    %         flag_QEC_TX_RX = 'TXLOQEC';
    %     end
    txwaveform = waveform;
    if size(txwaveform,1)==2
        % txwaveform = [I;Q]
        method_IMB_Demodulation = 'DirectIQDemod';
%         error(['Not support for method: ',method_IMB_Demodulation])
    elseif  size(txwaveform,1)==1
        % txwaveform = I-jQ
        method_IMB_Demodulation = 'DownConv2IQDemod';
    end
    
    % % generate LO for downconversion
    % [lo_I_Down_P,lo_Q_Down_P] = LO_WI_PN_IMB_g_0224(fc, fs, Nsamps, level_LO_dB, [], [], [], 'IQ', 0);
    
    %% 2020-03-16, Modify the filtering before DownConversion
    % FIR BPF or LPF
    FIR_f_cuttoff_tolerance = 0.01e6; % < fs/2
    FIR_Norder=10001;
    if size(bw_inband,2)==2
        FIR_Ftype='BPF';
        f_cutoff_L = bw_inband(1)+FIR_f_cuttoff_tolerance;
        f_cutof_H = bw_inband(2)-FIR_f_cuttoff_tolerance;
    elseif size(bw_inband,2)==1
        FIR_Ftype='LPF';
        f_cutoff_L = bw_inband(1)+FIR_f_cuttoff_tolerance;
        f_cutof_H = 0;
    end
    FIR_Wtype='REC';
    
    bw_fir_ratio_L = 2*f_cutoff_L/fs;
    FIR_w_cutoff_L = bw_fir_ratio_L*pi;
    bw_fir_ratio_H = 2*f_cutof_H/fs;
    FIR_w_cutoff_H = bw_fir_ratio_H*pi;
    
    % FIR generator
    b = fir_windowdesign_g(FIR_Norder, FIR_Ftype, FIR_w_cutoff_L, FIR_w_cutoff_H, FIR_Wtype);
    
    % FIR plot
    %     figure(2020022801)
    %     yyaxis right
    PLOT_FFT_dB_g(b*Nsamps, fs, Nsamps/1, [FIR_Ftype], 'df', 'full', 'pwr', [2020022801]);
    %     yyaxis left
    
    % apply FIR filtering
    %     txwaveform_FIR = conv(txwaveform,(b/2^0),'same'); % apply filter in time domain
    txwaveform_FIR = DSP_filter_g(b, txwaveform, 'FD');
    
    %     figure(2020022801)
    PLOT_FFT_dB_g(txwaveform_FIR, fs, Nsamps, ['wf wi LPF'], 1, 'full', 'pwr', [2020022801]);
    %     PLOT_FFT_dB_g(txwaveform, fs, Nsamps, ['wf wi LPF'], 1, 'full', 'pwr', [2020022801]);
    
    % downconversion for IQ seperation
    switch method_IMB_Demodulation
        case 'DirectIQDemod'
            wf_I = txwaveform_FIR(1,:);
            wf_Q = txwaveform_FIR(2,:);
            
        case 'DownConv2IQDemod'
            [wf_IQ_Down, wf_I_Down, wf_Q_Down] = Mixer_Up_Down_Convert_g(txwaveform_FIR, lo_I_IMB_Down, lo_Q_IMB_Down, fs, ['wf Downconversion for QEC'],'Down',[],[2020022801]);
            wf_I = wf_I_Down;
            wf_Q = wf_Q_Down;
            
%             wf_I = DSP_filter_g(b, wf_I_Down, 'FD');
%             wf_Q = DSP_filter_g(b, wf_Q_Down, 'FD');
    end
    
elseif strcmp(flag_QEC_TX_RX,'RXQEC')||isempty(lo_I_IMB_Down)||isempty(lo_Q_IMB_Down)
    %     flag_QEC_TX_RX = 'RXQEC';
    rxwaveform = waveform;
    wf_IQ = rxwaveform;
    wf_I = real(rxwaveform);
    wf_Q = imag(rxwaveform);
end

% plot
% figure(2020022801)
PLOT_FFT_dB_g(wf_I, fs, Nsamps, ['QEC wf I wi FIR'], 1, 'full', 'pwr', [2020022801]);
PLOT_FFT_dB_g(wf_Q, fs, Nsamps, ['QEC wf Q wi LPF'], 1, 'full', 'pwr', [2020022801]);
title('Spectrum for QEC application')

[IMB_MagDB, IMB_PhsDeg] = IQ_IMB_cor(wf_I, wf_Q,flag_QEC_TX_RX);
% if strcmp(method_IMB_Demodulation,'DirectIQDemod')
%     %     IMB_PhsDeg = -IMB_PhsDeg;
%     %     IMB_PhsDeg = IMB_PhsDeg-90; % V
%     if IMB_PhsDeg>0
%         IMB_PhsDeg = 90-IMB_PhsDeg; % V
%     elseif IMB_PhsDeg<0
%         IMB_PhsDeg = IMB_PhsDeg+90; % V
%     end
% elseif strcmp(method_IMB_Demodulation,'DownConv2IQDemod')
%     IMB_PhsDeg = IMB_PhsDeg; % compensation the DownConversion phase change
% end
IMB_MagDB_mean = mean(IMB_MagDB);
IMB_PhsDeg_mean = mean(IMB_PhsDeg);

%% 2020-03-12, Ignore !, Due to the filter has been applied.
% %% 2020-02-28, Tolearnce Error Detection
% if abs(IMB_MagDB)<=2 || abs(IMB_PhsDeg)<=2 % tolerance error
%     [IMB_MagDB_full, ~] = IQ_IMB_cor(wf_I, wf_Q,flag_QEC_TX_RX);
%     diff_IMB_MagDB =  IMB_MagDB_full-IMB_MagDB;
%     if abs(diff_IMB_MagDB)>=0.5
%         error('Tolearance Error! Please check the LO Phase Noise and SPURS')
%     else
%         IMB_MagDB = IMB_MagDB_full;
%     end
% end
%
% if abs(IMB_PhsDeg)<=2 % tolerance error
%     [~, IMB_PhsDeg_full] = IQ_IMB_cor(wf_I, wf_Q,flag_QEC_TX_RX);
%     diff_IMB_PhsDB =  IMB_PhsDeg_full-IMB_PhsDeg;
%     if abs(diff_IMB_PhsDB)>=0.5
%         error('Tolearance Error! Please check the LO Phase Noise and SPURS')
%     else
%         IMB_PhsDeg = IMB_PhsDeg_full;
%     end
% end

%% 2020-03-01, IQ waveform QEC correct
if strcmp(flag_QEC_TX_RX,'TXLOQEC')
    QEC_Mag = 10.^(0.5*IMB_MagDB_mean./20);
    QEC_Phs = 0.5*IMB_PhsDeg_mean./180*pi;
    I_Q_Matrix = [lo_I_IMB_Down/QEC_Mag(:);lo_Q_IMB_Down*QEC_Mag(:)]; % Magnitude correct
    IMB_Matrix = [cos(QEC_Phs(:)) sin(QEC_Phs(:)); sin(QEC_Phs(:)) +cos(QEC_Phs(:))]; % ???
    waveform_IQ_QEC(:,:) = IMB_Matrix\I_Q_Matrix; % phase correct
    
elseif strcmp(flag_QEC_TX_RX,'TXQEC') % TXQEC only support the estimation not correlation
    %     QEC_Mag = 10.^(0.5*IMB_MagDB./20);
    %     QEC_Phs = 0.5*IMB_PhsDeg./180*pi;
    %     for idBR=1:Nbr
    %         I_Q_Matrix = [txwaveform(idBR,:)/QEC_Mag(idBR);txwaveform(idBR,:)*QEC_Mag(idBR)]; % Magnitude correct
    %         IMB_Matrix = [cos(QEC_Phs(idBR)) sin(QEC_Phs(idBR)); -sin(QEC_Phs(idBR)) -cos(QEC_Phs(idBR))]; % IQ = I-1i*Q
    %         waveform_IQ_QECx = IMB_Matrix\I_Q_Matrix; % phase correct
    %         waveform_IQ_QEC(idBR,:) = complex(waveform_IQ_QECx(1,:),-waveform_IQ_QECx(2,:)); % IQ = I-1i*Q
    %     end
    waveform_IQ_QEC = [];
    if flag_debug_plot==1
        figure(20200514)
        PLOT_FFT_dB_g(waveform_IQ_QEC, fs, Nsamps, ['QEC wf'], 1, 'full');
    end
elseif strcmp(flag_QEC_TX_RX,'RXQEC')
    QEC_Mag = 10.^(0.5*IMB_MagDB./20);
    QEC_Phs = 0.5*IMB_PhsDeg./180*pi;
    for idBR=1:Nbr
        I_Q_Matrix = [wf_I(idBR,:)/QEC_Mag(idBR);wf_Q(idBR,:)*QEC_Mag(idBR)]; % Magnitude correct
        IMB_Matrix = [cos(QEC_Phs(idBR)) sin(QEC_Phs(idBR)); -sin(QEC_Phs(idBR)) -cos(QEC_Phs(idBR))]; % IQ = I-1i*Q
        waveform_IQ_QECx = IMB_Matrix\I_Q_Matrix; % phase correct
        waveform_IQ_QEC(idBR,:) = complex(waveform_IQ_QECx(1,:),-waveform_IQ_QECx(2,:)); % IQ = I-1i*Q
    end
else
    waveform_IQ_QEC = [];
end
% IMB_Matrix = [cos(QEC_Phs) -sin(QEC_Phs); -sin(QEC_Phs) +cos(QEC_Phs)]; % IQ = I+1i*Q
% IMB_Matrix = [cos(QEC_Phs) sin(QEC_Phs); -sin(QEC_Phs) -cos(QEC_Phs)]; % IQ = I-1i*Q

% export
if strcmp(flag_wf_original,'COLUMN')
    %     waveform_IQ_QEC=waveform_IQ_QEC(:); % COLUMN
    waveform_IQ_QEC=waveform_IQ_QEC.'; % COLUMN
end

if flag_debug_plot ==0
    %     close 2020080101
    close 2020022801
end
end

function [IMB_MagDB_est, IMB_PhsDeg_est] = IQ_IMB_cor(QECwf_I, QECwf_Q, flag_QEC_TX_RX)

% switch waveform to ROW
DIM_FFT = 2;

EyI = mean(QECwf_I.^2,DIM_FFT);
EyQ = mean(QECwf_Q.^2,DIM_FFT);
EyIyQ = mean(QECwf_I.*QECwf_Q,DIM_FFT);

IMB_MagDB_est = 10*log10(real(EyI./EyQ));
% flag_QEC_TX_RX = 'RXQEC';
switch flag_QEC_TX_RX
    case {'TXQEC'}
        IMB_PhsDeg_est = +asin(real(EyIyQ./(sqrt(EyI).*sqrt(EyQ))))/pi*180;
    case {'TXLOQEC'}
        IMB_PhsDeg_est = +asin(real(EyIyQ./(sqrt(EyI).*sqrt(EyQ))))/pi*180;
    case {'RXQEC'}
        IMB_PhsDeg_est = -asin(real(EyIyQ./(sqrt(EyI).*sqrt(EyQ))))/pi*180;
end
end