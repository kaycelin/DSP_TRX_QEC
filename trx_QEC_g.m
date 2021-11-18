%% 2020-02-28, Tolearnce Error Detection --> %% 2020-03-12, Ignore !, Due to the filter has been applied.
%% 2020-03-01, IQ waveform QEC correct
%% 2020-03-04, IQ waveform QEC correct for DL part ?
%% 2020-03-16, Modify the filtering before DownConversion
%% 2020-03-28, DSP_filter_g, Modify filtering form timedomain to FreqDomain
%% 2020-04-13, switch waveform to ROW: DIM_FFT = 2
%% 2021-03-11, txwaveform = [I_Row; Q_Row], and setup waveform's MAX_Size_Row=2
%% 2021-03-11, Multiple -1 to RX Image part



function [IMB_MagDB, IMB_PhsDeg, waveform_IQ_QEC] = trx_QEC_g3(waveform, lo_I_IMB_Down, lo_Q_IMB_Down, fs, bw_inband, flag_QEC_TX_RX, method_QECDemod, flag_debug_plot)
% flag_QEC_TX_RX = 'TXQEC'/'TXLOQEC'/'RXQEC'

Nsamps = length(waveform);
Nbr = 1;
df = fs/Nsamps;

% switch waveform to ROW
DIM_FFT = 2;
if size(waveform,1)>size(waveform,2) % COLUMN
    flag_wf_original = 'COLUMN';
    waveform=waveform.'; % switch to ROW
else
    flag_wf_original = 'ROW';
end

%% 2021-03-11, txwaveform = [I_Row; Q_Row], and setup waveform's MAX_Size_Row=2
if size(waveform,DIM_FFT-1)>2
    error('waveform MAX_Size_Row=2, [I_Row, Q_Row]!')
end

if isempty(lo_I_IMB_Down)
    lo_I_IMB_Down = 1;
end

if isempty(lo_Q_IMB_Down)
    lo_Q_IMB_Down = 1;
end

if any([lo_I_IMB_Down,lo_Q_IMB_Down])==1
    disp_DnConv_QEC = ['wf for QEC'];
else
    disp_DnConv_QEC = ['wf Downconversion for QEC'];
end

if size(lo_I_IMB_Down,1)>size(lo_I_IMB_Down,2) % COLUMN
    lo_I_IMB_Down=lo_I_IMB_Down.'; % switch to ROW
    lo_Q_IMB_Down=lo_Q_IMB_Down.'; % switch to ROW
end

if ~exist('flag_debug_plot','var')||isempty(flag_debug_plot)
    flag_debug_plot = 0;
end

%% 2020-03-16, Modify the filtering before DownConversion
if exist('bw_inband','var')&&~isempty(bw_inband)
    % FIR BPF or LPF
    FIR_f_cuttoff_tolerance = 0.01e6; % < fs/2
    FIR_Norder=10001;
    if size(bw_inband,2)==1
        FIR_Ftype='LPF';
        f_cutoff_L = fix((bw_inband(1) + sign(bw_inband(1))*FIR_f_cuttoff_tolerance)/df)*df;
        f_cutoff_H = 0;
    elseif size(bw_inband,2)==2
        FIR_Ftype='BPF';
        f_cutoff_L = fix((bw_inband(1) - sign(bw_inband(1))*FIR_f_cuttoff_tolerance)/df)*df;
        if f_cutoff_L<0
            f_cutoff_L=0;
        end
        f_cutoff_H = fix((bw_inband(2) + sign(bw_inband(2))*FIR_f_cuttoff_tolerance)/df)*df;
    else
        error('!')
    end
    
    FIR_Wtype='REC';
    
    bw_fir_ratio_L = 2*f_cutoff_L/fs;
    FIR_w_cutoff_L = bw_fir_ratio_L*pi;
    bw_fir_ratio_H = 2*f_cutoff_H/fs;
    FIR_w_cutoff_H = bw_fir_ratio_H*pi;
    
    % FIR generator
    b = fir_windowdesign_g(FIR_Norder, FIR_Ftype, FIR_w_cutoff_L, FIR_w_cutoff_H, FIR_Wtype);
    
    % FIR plot
    PLOT_FFT_dB_g(b*Nsamps, fs, Nsamps/1, [FIR_Ftype], 'df', 'full', 'pwr', [2020022801]);
    
else
    b = [];
end

if ~exist('flag_QEC_TX_RX','var')||isempty(flag_QEC_TX_RX)||strcmp(flag_QEC_TX_RX,'TXQEC')||strcmp(flag_QEC_TX_RX,'TXLOQEC')
    txwaveform = waveform;
    
    if ~exist('method_QECDemod','var')||isempty(method_QECDemod)
        error('method_QECDemod?')
    elseif strcmp(method_QECDemod,'Direct_IQDemod')
        method_FIR = 'on'; % default
        if strcmp(method_FIR,'on')&&~isempty(b)
            txwaveform_FIR = DSP_filter_g(b, txwaveform, 'FD');
            disp_FIR = ' wi BPF';
        else
            txwaveform_FIR = txwaveform;
            disp_FIR = [];
        end
        
        txwaveform_FIR = DSP_filter_g(b, txwaveform, 'FD');
        wf_I = txwaveform_FIR(1,:);
        wf_Q = txwaveform_FIR(2,:);
        
        [IMB_MagDB, IMB_PhsDeg] = IQ_IMB_cor(wf_I, wf_Q,flag_QEC_TX_RX)
        
    elseif strcmp(method_QECDemod,'DNConv_FIR_IQDemod')
        % downconversion for IQ seperation
        method_ExchangeLOIQ = 'off';
        %         if strcmp(method_ExchangeLOIQ,'off')
        if size(txwaveform,DIM_FFT-1)==2
            [wf_IQ_Down, wf_I_Down, wf_Q_Down] = Mixer_Up_Down_Convert_g(txwaveform, lo_I_IMB_Down, lo_Q_IMB_Down, fs, [disp_DnConv_QEC],'Down',[],[2020022801]);
            %         elseif strcmp(method_ExchangeLOIQ,'on')
            %             [wf_IQ_Down, wf_I_Down, wf_Q_Down] = Mixer_Up_Down_Convert_g(txwaveform, lo_Q_IMB_Down, lo_I_IMB_Down, fs, ['wf Downconversion for QEC'],'Down',[],[2020022801]);
        else
            wf_I_Down = real(txwaveform);
            wf_Q_Down = -1*imag(txwaveform);
            %% 2021-03-11, Multiple -1 to RX Image part
        end

        method_FIR = 'on'; % default
        if strcmp(method_FIR,'on')&&~isempty(b)
            wf_I = DSP_filter_g(b, wf_I_Down, 'FD');
            wf_Q = DSP_filter_g(b, wf_Q_Down, 'FD');
            disp_FIR = ' wi LPF';
        else
            wf_I = wf_I_Down;
            wf_Q = wf_Q_Down;
            disp_FIR = [];
        end
        %         PLOT_FFT_dB_g(b*Nsamps, fs, Nsamps/1, [FIR_Ftype], 'df', 'full', 'pwr', [2020022801]);
        %         PLOT_FFT_dB_g(wf_I, fs, Nsamps, ['QEC wf I wi FIR'], 1, 'full', 'pwr', [2020022801]);
        
        if strcmp(method_ExchangeLOIQ,'off')
            [IMB_MagDB, IMB_PhsDeg] = IQ_IMB_cor(wf_I, wf_Q,flag_QEC_TX_RX)
            %         elseif strcmp(method_ExchangeLOIQ,'on')
            %             [IMB_MagDB, IMB_PhsDeg] = IQ_IMB_cor(wf_I, wf_Q,flag_QEC_TX_RX)
        end
        
    elseif strcmp(method_QECDemod,'DNConvLOI_FIR_IQDemod')
        if any([~isempty(lo_I_IMB_Down),~isempty(lo_Q_IMB_Down)])
            [wf_IQ_Down, wf_I_Down, wf_Q_Down] = Mixer_Up_Down_Convert_g(txwaveform, lo_I_IMB_Down, lo_I_IMB_Down, fs, [disp_DnConv_QEC],'Down',[],[2020022801]);
        else
            error('?')
        end
        
        method_FIR = 'on'; % default
        if strcmp(method_FIR,'on')&&~isempty(b)
            wf_I = DSP_filter_g(b, wf_I_Down, 'FD');
            wf_Q = DSP_filter_g(b, wf_Q_Down, 'FD');
            disp_FIR = ' wi LPF';
        else
            wf_I = wf_I_Down;
            wf_Q = wf_Q_Down;
            disp_FIR = [];
        end
        
        [IMB_MagDB, IMB_PhsDeg] = IQ_IMB_cor(wf_I, wf_Q,flag_QEC_TX_RX)
    end
    
elseif strcmp(flag_QEC_TX_RX,'RXQEC')||isempty(lo_I_IMB_Down)||isempty(lo_Q_IMB_Down)
    rxwaveform = waveform;
    if size(waveform,DIM_FFT-1)==2 % [I_Row, Q_Row]
        wf_I_RX = rxwaveform(1,:);
        wf_Q_RX = -1*rxwaveform(2,:);
        %% 2021-03-11, Multiple -1 to RX Image part
    else % [I+jQ]
        wf_I_RX = real(rxwaveform);
        wf_Q_RX = imag(rxwaveform); 
    end
    
    method_FIR = 'on'; % default
    if strcmp(method_FIR,'on')&&~isempty(b)
        wf_I = DSP_filter_g(b, wf_I_RX, 'FD');
        wf_Q = DSP_filter_g(b, wf_Q_RX, 'FD');
        disp_FIR = ' wi LPF';
    else
        wf_I = wf_I_RX;
        wf_Q = wf_Q_RX;
        disp_FIR = [];
    end
   
    [IMB_MagDB, IMB_PhsDeg] = IQ_IMB_cor(wf_I, wf_Q,flag_QEC_TX_RX)
    
end

PLOT_FFT_dB_g(wf_I, fs, Nsamps, ['QEC wf I',disp_FIR], 1, 'full', 'pwr', [2020022801]);
PLOT_FFT_dB_g(wf_Q, fs, Nsamps, ['QEC wf Q',disp_FIR], 1, 'full', 'pwr', [2020022801]);
title(['Spectrum of ', flag_QEC_TX_RX,' Analysis, IMBMagDB=',num2str(round(IMB_MagDB,2)),' IMBPhsDeg=',num2str(round(IMB_PhsDeg,2))])

% [IMB_MagDB, IMB_PhsDeg] = IQ_IMB_cor(wf_I, wf_Q,flag_QEC_TX_RX);
% IMB_MagDB_mean = mean(IMB_MagDB);
% IMB_PhsDeg_mean = mean(IMB_PhsDeg);

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
    %     QEC_Mag = 10.^(0.5*IMB_MagDB./20);
    %     QEC_Phs = 0.5*IMB_PhsDeg./180*pi;
    %     I_Q_Matrix = [lo_I_IMB_Down/QEC_Mag(:);lo_Q_IMB_Down*QEC_Mag(:)]; % Magnitude correct
    %     IMB_Matrix = [cos(QEC_Phs(:)) sin(QEC_Phs(:)); sin(QEC_Phs(:)) +cos(QEC_Phs(:))]; % ???
    %     waveform_IQ_QEC(:,:) = IMB_Matrix\I_Q_Matrix; % phase correct
    
    QEC_Mag = 10.^(0.5*IMB_MagDB./20);
    QEC_Phs = 0.5*IMB_PhsDeg./180*pi;
    I_Q_Matrix = [txwaveform(1,:)/QEC_Mag(:);txwaveform(2,:)*QEC_Mag(:)]; % Magnitude correct
    IMB_Matrix = [cos(QEC_Phs(:)) sin(QEC_Phs(:)); sin(QEC_Phs(:)) +cos(QEC_Phs(:))]; % ???
    IMB_Matrix = [cos(QEC_Phs(:)) sin(QEC_Phs(:)); -sin(QEC_Phs(:)) -cos(QEC_Phs(:))]; % IQ = I-1i*Q
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
%     if flag_debug_plot==1
%         PLOT_FFT_dB_g(waveform_IQ_QEC, fs, Nsamps, ['QEC wf'], 1, 'full','pwr',20200514);
%     end
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
    close 2020022801
end
end

function [IMB_MagDB_est, IMB_PhsDeg_est] = IQ_IMB_cor(QECwf_I, QECwf_Q, flag_QEC_TX_RX)

% switch waveform to ROW
DIM_FFT = 2;

EyI = mean(QECwf_I.^2,DIM_FFT);
EyQ = mean(QECwf_Q.^2,DIM_FFT);
EyIyQ = mean(QECwf_I.*QECwf_Q,DIM_FFT);
IMB_MagDB_est = 20*log10(sqrt(EyI./EyQ));

% EyIyQ2 = (QECwf_I.*QECwf_Q);
% EyI2 = (QECwf_I.^2);
% EyQ2 = (QECwf_Q.^2);
% flag_QEC_TX_RX = 'RXQEC';
switch flag_QEC_TX_RX
    case {'TXQEC'}
        IMB_PhsDeg_est = +asin(real(EyIyQ./(sqrt(EyI).*sqrt(EyQ))))/pi*180;
        %         IMB_PhsDeg_est = +asin((real(EyIyQ2./((sqrt(EyI2+EyQ2))))))/pi*180;
    case {'TXLOQEC'}
        IMB_PhsDeg_est = +asin(real(EyIyQ./(sqrt(EyI).*sqrt(EyQ))))/pi*180;
    case {'RXQEC'}
        IMB_PhsDeg_est = -asin(real(EyIyQ./(sqrt(EyI).*sqrt(EyQ))))/pi*180;
end
end