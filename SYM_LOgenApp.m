function [lo_ideal,lo_realistic,tableInputLO,tableOutputLO] = LOgenApp(LO, fs, Nsamps, flag_LO, fnum, flag_LO_plot, disp_title, fnum_save_dir)
%% 2020-09-15, LO: -leveldB/-fLO/-PN/-IMB/-SPURS
%% 2020-09-15, flag_LO: -PhsNoise/-IMB/-SPURS/-AMtoPMPhsDriftDeg/-QEC
%% 2020-09-15, method_QEC: 'DNConv_FIR_IQDemod'/'Direct_IQDemod'
%% 2020-09-16, QEC demodulation signal should NOT based on ZERO_HZ
%% 2021-02-28, Update

df = fs/Nsamps;

if ~exist('flag_LO_plot','var')||isempty(flag_LO_plot)
    flag_LO_plot = 'semilogx';
elseif strcmp(flag_LO_plot,'semilogx')||strcmp(flag_LO_plot,'full')||strcmp(flag_LO_plot,'half')
    flag_LO_plot = flag_LO_plot;
else
    error('!')
end

if ~exist('disp_title','var')||isempty(disp_title)
    disp_title=[];
else
    disp_title = [', ',disp_title];
end


%% LO Generator:
%% LO Generated Perfactly!!
if isstruct(LO)
    if isfield(LO,'leveldB')
        LO_leveldB = LO.leveldB;
    else
        LO_leveldB = 0;
    end
    if isfield(LO,'fLO')
        fLO = ceil(LO.fLO/df)*df;
    else
        error('Input of fLO!')
    end
    
    [lo_I_ideal, lo_Q_ideal] = LO_WI_PN_IMB_g(fLO, fs, Nsamps, LO_leveldB, [], [], [], 'IQ', [],0);
    lo_ideal = [lo_I_ideal;lo_Q_ideal];
else
    LO_leveldB = 0;
    fLO = LO;
end

disp(['Generate LO with ',num2str(fLO/1e6),'MHz ===================================================='])
tableInputLO_LevelandFreq = table(LO_leveldB, fLO)
[lo_I_ideal, lo_Q_ideal] = LO_WI_PN_IMB_g(fLO, fs, Nsamps, LO_leveldB, [], [], [], 'IQ', [],0);
lo_ideal = [lo_I_ideal;lo_Q_ideal];

if exist('fnum','var')&&~isempty(fnum)
    %     figure(fnum*200)
    %     subplot(1,2,1)
    LO_I=PLOT_FFT_dB_g(lo_I_ideal, fs, Nsamps, ['LO I ideal'], 'Hz', flag_LO_plot, 'pwr', [200*fnum(1),1,2,1]);
    %     LO_I=PLOT_FFT_dB_g(lo_I_Up, fs, Nsamps, ['LO Up'], 'Hz', 'half');
    ylim([-400, LO_leveldB+10]), title(['LO I',disp_title])
    
    %     subplot(1,2,2)
    LO_Q=PLOT_FFT_dB_g(lo_Q_ideal, fs, Nsamps, ['LO Q ideal'], 'Hz', flag_LO_plot, 'pwr', [200*fnum(1),1,2,2]);
    ylim([-400, LO_leveldB+10]), title(['LO Q',disp_title])
end

% pwr check
[IpwrdB_lo_I_ideal] = Pwr_Inband_g(fft(lo_I_ideal), fs, [fLO fLO], 0, 'half', 0);

if ~exist('flag_LO','var')||isempty(flag_LO)
    lo_realistic = [];
    tableInputLO = [];
    tableOutputLO = [];
    return
else
    %% LO Generated with Phase Noise/IQ Imbalance/SPURS
    if isfield(flag_LO,'PhsNoise')
        flag_LO_PhsNoise = flag_LO.PhsNoise;
    else
        flag_LO_PhsNoise = 0;
    end
    if isfield(flag_LO,'IMB')
        flag_LO_IMB = flag_LO.IMB;
    else
        flag_LO_IMB = 0;
    end
    if isfield(flag_LO,'SPURS')
        flag_LO_SPURS = flag_LO.SPURS;
    else
        flag_LO_SPURS = 0;
    end
    if isfield(flag_LO,'AMtoPMPhsDriftDeg')
        flag_LO_AMtoPMPhsDriftDeg = flag_LO.AMtoPMPhsDriftDeg;
    else
        flag_LO_AMtoPMPhsDriftDeg = 0;
    end
    
    %     if ~isempty(flag_LO.AMtoPMPhsDriftDeg)&&(flag_LO.AMtoPMPhsDriftDeg~=0)
    %         flag_LO_AMtoPMPhsDriftDeg = flag_LO.AMtoPMPhsDriftDeg;
    %     else
    %         flag_LO_AMtoPMPhsDriftDeg = 0;
    %     end
    % if isfield(flag_LO,'QEC')&&( (ischar(flag_LO.QEC)&&~strcmp(flag_LO.QEC,'off'))||(isscalar(flag_LO.QEC)&&(flag_LO.QEC~=0)) )
    if isfield(flag_LO,'QEC')
        flag_LO_QEC = flag_LO.QEC;
    else
        flag_LO_QEC = 'off';
    end
    tableInputLO_flag  = table(flag_LO_PhsNoise, flag_LO_IMB, flag_LO_SPURS , flag_LO_AMtoPMPhsDriftDeg, flag_LO_QEC);
end

if isfield(LO,'PN')
    if flag_LO_PhsNoise == 0
        LO.PN = [];
        disp_LO_PN = ['LO wo Phase Noise'];
        tableInputLO_PN = 'NA';
    else
        Length_f_offset_Hz = length(LO.PN.f_offset_Hz);
        Length_g_offset_dBc1Hz = length(LO.PN.g_offset_dBc1Hz);
        Max_Length_PN = max(Length_f_offset_Hz,Length_g_offset_dBc1Hz)
        
        LO.PN.f_offset_Hz(end) = fs/2-0; % phase noise spectrum, frequencies
        %% 2021-03-02, add eps ???
        PN_f_offset_Hz = repmat(LO.PN.f_offset_Hz(end),[Max_Length_PN,1])+repmat(eps, [Max_Length_PN,1]).*10.^(1:Max_Length_PN)';
        PN_g_offset_dBc1Hz = repmat(LO.PN.g_offset_dBc1Hz(end),[Max_Length_PN,1])+repmat(eps, [Max_Length_PN,1]).*10.^(1:Max_Length_PN)';
        
        PN_f_offset_Hz(1:Length_f_offset_Hz) = LO.PN.f_offset_Hz(:);
        PN_g_offset_dBc1Hz(1:Length_g_offset_dBc1Hz) = LO.PN.g_offset_dBc1Hz(:);
        
        tableInputLO_PN = {table(PN_f_offset_Hz, PN_g_offset_dBc1Hz)};
        
        LO.PN.f_offset_Hz = PN_f_offset_Hz;
        LO.PN.g_offset_dBc1Hz = PN_g_offset_dBc1Hz;
        
        if flag_LO_PhsNoise == 1
            disp_LO_PN = ['LO wi Phase Noise'];
        elseif flag_LO_PhsNoise > 1
            disp_LO_PN = ['LO wi Phase Noise, delta:',num2str(flag_LO_PhsNoise),'dB'];
        end
    end
else
    LO.PN = [];
    disp_LO_PN = ['LO wo Phase Noise'];
    tableInputLO_PN = 'NA';
end

% Mixer IQ Imbalance setup: IMB_Phs_deg ; IMB_Mag_dB
if isfield(LO,'IMB')
    if flag_LO_IMB == 0
        LO.IMB = [];
        disp_LO_IMB = ['LO wo IQ imbalance'];
        tableInputLO_IMB = 'NA';
    elseif flag_LO_IMB == 1
        IMB_PhsDeg = LO.IMB.PhsDeg;
        IMB_MagdB = LO.IMB.MagdB;
        tableInputLO_IMB =  table(IMB_MagdB, IMB_PhsDeg);
        disp_LO_IMB = ['Phs IMB:', num2str(IMB_PhsDeg), 'deg', '; Mag IMB:', num2str(IMB_MagdB), 'dB'];
    end
else
    LO.IMB = [];
    disp_LO_IMB = ['LO wo IQ imbalance'];
    tableInputLO_IMB = 'NA';
end

if isfield(LO,'SPURS')
    if flag_LO_SPURS == 0
        LO.SPURS = [];
        disp_LO_SPURS = ['LO wo SPURS'];
        tableInputLO_SPURS = 'NA';
    elseif flag_LO_SPURS == 1
        disp_LO_SPURS = ['LO wi SPURS'];
        foffset_spurs_Hz = LO.SPURS.foffset_spurs_Hz(:);
        g_spurs_dBc1Hz = LO.SPURS.g_spurs_dBc1Hz(:);
        tableInputLO_SPURS =  {table(foffset_spurs_Hz, g_spurs_dBc1Hz)};
    elseif flag_LO_SPURS > 1
        foffset_spurs_Hz = LO.SPURS.foffset_spurs_Hz(:);
        g_spurs_dBc1Hz = LO.SPURS.g_spurs_dBc1Hz(:);
        %     disp_LO_SPURS = ['LO wi SPURS, delta:',num2str(flag_LO_SPURS),'dB'];
        disp_LO_SPURS = ['LO wi SPUR*',num2str(length(g_spurs_dBc1Hz)),', MaxSPUR:',num2str(max(g_spurs_dBc1Hz)),'dBc'];
        tableInputLO_SPURS =  {table(foffset_spurs_Hz, g_spurs_dBc1Hz)};
    end
else
    LO.SPURS = [];
    disp_LO_SPURS = ['LO wo SPURS'];
    tableInputLO_SPURS = 'NA';
end

% export table
tableInputLO = table(tableInputLO_LevelandFreq, tableInputLO_flag, tableInputLO_PN, tableInputLO_IMB, tableInputLO_SPURS)


%% Generate LO with Phase Noise and IQ Imbalance
[lo_I_PN_IMB_SPURS,lo_Q_PN_IMB_SPURS] = LO_WI_PN_IMB_g(fLO, fs, Nsamps, LO.leveldB, LO.PN, LO.IMB, LO.SPURS, 'IQ', flag_LO_AMtoPMPhsDriftDeg,0);
% [lo_I_Up_PN_IMB_SPURS,lo_Q_Up_PN_IMB_SPURS] = LO_WI_PN_IMB_g3(LOU.fLO, fs, Nsamps, LOU.leveldB, [], IMB_LOU, [], 'IQ', flag_LOU.AMtoPMPhsDriftDeg,0);

% plot
disp_LO_PN_IMB_SPURS_Up = [disp_LO_PN, sprintf('\n'), disp_LO_IMB, sprintf('\n'), disp_LO_SPURS];

if exist('fnum','var')&&~isempty(fnum)
    %     figure(fnum*200)
    %     subplot(1,2,1)
    PLOT_FFT_dB_g((lo_I_PN_IMB_SPURS), fs, length(lo_I_PN_IMB_SPURS), [disp_LO_PN_IMB_SPURS_Up], 'Hz', flag_LO_plot, 'pwr', [200*fnum(1),1,2,1]);
    title(['LO I',disp_title])
    
    %     subplot(1,2,2)
    PLOT_FFT_dB_g(lo_Q_PN_IMB_SPURS, fs, length(lo_Q_PN_IMB_SPURS), [disp_LO_PN_IMB_SPURS_Up], 'Hz', flag_LO_plot, 'pwr', [200*fnum(1),1,2,2]);
    title(['LO Q',disp_title])
end

%% Check LO Phase Noise Pwr
[tableOutputLO_PN{1,1}, tableOutputLO_SPURS{1,1}, tableOutputLO_IMB] = LO_Performance_g(lo_I_PN_IMB_SPURS, lo_Q_PN_IMB_SPURS, fs, LO);

if isempty(tableOutputLO_PN)
    disp(disp_LO_PN)
    tableOutputLO_PN='NA';
else
    tableOutputLO_PN{:}
end
if isempty(tableOutputLO_SPURS)
    disp(disp_LO_SPURS)
    tableOutputLO_SPURS='NA';
else
    tableOutputLO_SPURS{:}
end
if isempty(tableOutputLO_IMB)
    disp(disp_LO_IMB)
    tableOutputLO_IMB='NA';
else
    tableOutputLO_IMB
end

%% QEC for the Imbalance LO
% if ( ischar(flag_LO_QEC)&&strcmp(flag_LO_QEC,'on')||(isscalar(flag_LO_QEC)&&flag_LO_QEC~=0) )&&(flag_LO_IMB==1) % 1:ON; 0:OFF
if ~any(flag_LO_QEC==0)&&~strcmp(flag_LO_QEC,'off')
    %% 2020-09-15, method_QEC: 'DNConv_FIR_IQDemod'/'Direct_IQDemod'
    %% 2020-09-16, QEC demodulation signal should NOT based on ZERO_HZ
    fcuttoff_QEC = 5e6;
    if ~strcmp(flag_LO_QEC,'RXLOQEC')
        method_QEC = 'DNConv_FIR_IQDemod';
        %         method_QEC = 'DNConvLOI_FIR_IQDemod';
        %         method_QEC = 'DNConv_FIRDirect_IQDemod_IQDemod';
        switch method_QEC
            case 'DNConv_FIR_IQDemod'
                %% 1. Generate QEC CW source from fNCO=10MHz
                fBBQEC = 10e6;
                t = (0:Nsamps-1)/fs;
                TxLoQECBBCSource = exp(1i*2*pi*fBBQEC*t);
                
                %% 2. Mix with fLOU
                %% 2.1. Generate txQECSource by Mix with Imbalance fLOU
                %% 2.2. Filter the Image part by FIR[bwRFQEC]
                bwRFQEC = fix((fBBQEC+fLO+fcuttoff_QEC*[-1 1])/df)*df;
                method_QEC_MixerUp = 'wiFIR';
                if strcmp(method_QEC_MixerUp, 'wiFIR')
                    [TxLoQECRFSource, TxLoQECRFSourceI, TxLoQECRFSourceQ] = Mixer_Up_Down_Convert_g(TxLoQECBBCSource, lo_I_PN_IMB_SPURS, lo_Q_PN_IMB_SPURS, fs, ['txQECSource'],'Up',[bwRFQEC],[0803]);
                elseif strcmp(method_QEC_MixerUp, 'woFIR')
                    error(['method_QEC_MixerUp: ',method_QEC_MixerUp])
                    [TxLoQECRFSource, TxLoQECRFSourceI, TxLoQECRFSourceQ] = Mixer_Up_Down_Convert_g(TxLoQECBBCSource, lo_I_PN_IMB_SPURS, lo_Q_PN_IMB_SPURS, fs, ['txQECSource'],'Up',[],[0803]);
                end
                
                %% 3. QEC from fBB_UpConversionSource
                %% 3.1. Downconvert to fBB with fLOU by ideal LO
                %% 3.2. Filter the HighFreqBand by FIR[bwBBQEC]
                % fBB, DNConv_FIR_IQDemod
                bwBBQEC = fix((0+0+(fBBQEC+fcuttoff_QEC)*[-1 1])/df)*df;
                [QECest_lo_MagdB, QECest_lo_PhsDeg, ~] = trx_QEC_g([TxLoQECRFSourceI;TxLoQECRFSourceQ], lo_I_ideal, lo_Q_ideal, fs, bwBBQEC,'TXLOQEC', 'DNConv_FIR_IQDemod');
                
            case 'DNConvLOI_FIR_IQDemod'
                %% 1. QEC from Original Realistic LO Source
                %% 1.1. Convert to 2*fLO with fLOU by ideal LO
                %% 1.2. Filter the LowFreqBand by FIR[bwBBQEC]
                % fBBQEC, DNConv_FIR_IQDemod
                fBBQEC = 10e6;
                fLODQEC=fLO-fBBQEC;
                t = [1:Nsamps]/fs;
                bwBBQEC = fix((0+0+(fBBQEC+fcuttoff_QEC)*[-1 1])/df)*df;
                %             [QECest_lo_MagdB, QECest_lo_PhsDeg, loQECwf] = trx_QEC_g([lo_I_PN_IMB_SPURS;lo_Q_PN_IMB_SPURS], sin(2*pi*fLO_DN*t), sin(2*pi*fLO_DN*t), fs, bwBBQEC,'TXLOQEC', 'DNConv_LOI_FIR_IQDemod');
                [QECest_lo_MagdB, QECest_lo_PhsDeg, ~] = trx_QEC_g([lo_I_PN_IMB_SPURS;lo_Q_PN_IMB_SPURS], cos(2*pi*fLODQEC*t), cos(2*pi*fLODQEC*t), fs, bwBBQEC,'TXLOQEC', 'DNConvLOI_FIR_IQDemod');
                
            case 'Direct_IQDemod'
                %% 1. QEC from Original Realistic LO Source
                %% 1.1. Demodulation directly, suitable for RFSampling RX
                %% 1.2. No Filter Needed
                % fLO, Direct_IQDemod
                bwLOQEC = fix((fLO+0+fcuttoff_QEC*[-1 1])/df)*df;
                %             [QECest_lo_MagdB, QECest_lo_PhsDeg, loQECwf] = trx_QEC_g([lo_I_PN_IMB_SPURS;lo_Q_PN_IMB_SPURS], [], [], fs, bwLOQEC,'TXLOQEC', 'Direct_IQDemod');
                [QECest_lo_MagdB, QECest_lo_PhsDeg, ~] = trx_QEC_g([lo_I_PN_IMB_SPURS;lo_Q_PN_IMB_SPURS], [], [], fs, [],'TXLOQEC', 'Direct_IQDemod');
        end
    elseif strcmp(flag_LO_QEC,'RXLOQEC')
        bwRFQEC = fix((fBBQEC+fLO+fcuttoff_QEC*[-1 1])/df)*df;
        %% Generate RxLoQECRFSourceI by QEC Mixer for flag_LO_QEC=RXLOQEC
        lo_QEC_RXLOQEC = cos(2*pi*fBBQEC*t);
        [~, RxLoQECRFSourceI, ~] = Mixer_Up_Down_Convert_g(lo_I_PN_IMB_SPURS, lo_QEC_RXLOQEC, lo_QEC_RXLOQEC, fs, ['rxQECSource'],'Up',[bwRFQEC],[0803]);
        %% 2021-03-09, Question1, empty Q path for lo_QEC_RXLOQEC? ...
        %% [~, RxLoQECRFSourceI, ~] = Mixer_Up_Down_Convert_g(lo_I_PN_IMB_SPURS, lo_QEC_RXLOQEC, [], fs, ['rxQECSource'],'Up',[bwRFQEC],[0803]);
        %% 2021-03-09, Question2, Get lo_QEC_RXLOQEC for TX Mixer I path ?
    
        method_FIRBBQEC = 'on';
        switch method_FIRBBQEC
            case 'on'
                bwBBQEC = fix((0+0+(fBBQEC+fcuttoff_QEC)*[-1 1])/df)*df;
            case 'off'
                bwBBQEC = [];
        end
        [RxLoBBQECSource, RxLoBBQECSourceI, RxLoBBQECSourceQ] = Mixer_Up_Down_Convert_g(RxLoQECRFSourceI, lo_I_PN_IMB_SPURS, lo_Q_PN_IMB_SPURS, fs, ['BBQECSource'],'Down',[bwBBQEC],[0803]);
        [QECest_MagdB, QECest_PhsDeg, ~] = trx_QEC_g(RxLoBBQECSource, [], [], fs, [bwBBQEC], 'RXQEC','DNConv_FIR_IQDemod');
        
    end
    close 0803
    %     figure(20200915), plot(lo_I_PN_IMB_SPURS), hold on, plot(lo_I_ideal)
    tableOutput_LOQEC_Estimation = table(QECest_lo_MagdB, QECest_lo_PhsDeg);
    
    % generate the LO with QEC(Quadurature Error Correction)
    %% 2020-09-17, TXLO QEC
    QEC_Mag = 10.^(0.5*QECest_lo_MagdB./20);
    QEC_Phs = 0.5*QECest_lo_PhsDeg./180*pi;
    I_Q_Matrix = [lo_I_PN_IMB_SPURS/QEC_Mag(:);lo_Q_PN_IMB_SPURS*QEC_Mag(:)]; % Magnitude correct
    IMB_Matrix = [cos(QEC_Phs(:)) sin(QEC_Phs(:)); sin(QEC_Phs(:)) +cos(QEC_Phs(:))]; % ???
    waveform_IQ_QEC(:,:) = IMB_Matrix\I_Q_Matrix; % phase correct
    
    loQEC_I = waveform_IQ_QEC(1,:);
    loQEC_Q = waveform_IQ_QEC(2,:);
    
    %% 4. Verify the New txwf by loQEC
    switch method_QEC
        case 'Direct_IQDemod'
            [QECest_loQEC_MagdB, QECest_loQEC_PhsDeg, ~] = trx_QEC_g([loQEC_I;loQEC_Q], [], [], fs, [],'TXLOQEC', 'Direct_IQDemod');
        case 'DNConv_LOI_FIR_IQDemod'
            [QECest_loQEC_MagdB, QECest_loQEC_PhsDeg, ~] = trx_QEC_g([loQEC_I;loQEC_Q], cos(2*pi*fLODQEC*t), cos(2*pi*fLODQEC*t), fs, bwBBQEC,'TXLOQEC', 'DNConv_LOI_FIR_IQDemod');
        case 'DNConv_FIR_IQDemod'
            [txQECwfVerify, ~, ~] = Mixer_Up_Down_Convert_g(TxLoQECRFSource, loQEC_I, loQEC_Q, fs, ['txQECwfVerify'],'Down',[bwBBQEC],[0803]);
            [QECest_loQEC_MagdB, QECest_loQEC_PhsDeg, ~] = trx_QEC_g(txQECwfVerify,  loQEC_I, loQEC_Q, fs, bwBBQEC,'RXQEC');
    end
    close 0803
    tableOutput_LOQEC_Correction = table(QECest_loQEC_MagdB, QECest_loQEC_PhsDeg);
    
    % export LO
    disp_lo_QEC = ['QEClo'];
    lo_I_Realistic = real(loQEC_I);
    lo_Q_Realistic = real(loQEC_Q);
    tableOutputLO = table(tableOutputLO_PN, tableOutputLO_SPURS, tableOutputLO_IMB, tableOutput_LOQEC_Estimation, tableOutput_LOQEC_Correction)
    
    % plot
    disp_QEC = ['LO with QEC'];
    PLOT_FFT_dB_g((lo_I_Realistic), fs, length(lo_I_Realistic), [disp_QEC], 'Hz', flag_LO_plot, 'pwr', [200*fnum(1),1,2,1]);
    title(['LO I',disp_title])
    PLOT_FFT_dB_g((lo_Q_Realistic), fs, length(lo_Q_Realistic), [disp_QEC], 'Hz', flag_LO_plot, 'pwr', [200*fnum(1),1,2,2]);
    title(['LO Q',disp_title])
    
else
    disp_lo_QEC = [];
    lo_I_Realistic = real(lo_I_PN_IMB_SPURS);
    lo_Q_Realistic = real(lo_Q_PN_IMB_SPURS);
    tableOutputLO = table(tableOutputLO_PN, tableOutputLO_SPURS, tableOutputLO_IMB)
end
lo_realistic=[lo_I_Realistic;lo_Q_Realistic];

if exist('fnum','var')&&~isempty(fnum)
    %% 2020-10-31, fnum_save_dir: save picture to folder
    if exist('fnum_save_dir')&&~isempty(fnum_save_dir)
        fnum_save_file = [fnum_save_dir,'\',num2str(fnum),'.fig']
        saveas(gcf,[fnum_save_file])
    end
end

end