%% 2021-03-10, PA1

%% 1. RX Mixer QEC *******************************************************
%% 1a. Generate real signal of RXQECBBCSource_real, 10MHz
% Input 1a: ========================================
fs = 3.9322e+09
Nsamps = 2457600
df = fs/Nsamps
fBBQEC = 10e6
fnum = 111801
% ========================================

t = (0:Nsamps-1)/fs;
RXQECBBSource_real = real(exp(1i*2*pi*fBBQEC*t));
PLOT_FFT_dB_g(RXQECBBSource_real, fs, Nsamps, 'RXQECBBSource Real', 'df', 'full', 'pwr', fnum, [], []);
fnum = fnum+1

%% 1b. Generate LOU for Real Signal UpConversion
% Input 1b: ========================================
edit SystemSim_LoadCoefficient.m

xlsFile='Matlab_Linkbudget.xlsx';
xlsSheetName = 'QEC';
xlsRangeArrayInput ='A1:B550';
xlsReadShift = [0 0];
xlsWriteShift = [0 0];

LOU_AMtoPMPhsDriftDeg = ExcelRead('LOU_AMtoPMPhsDriftDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOU_leveldB = ExcelRead('LOU_leveldB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOU_fLO = ExcelRead('LOU_fLO',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOU_PN_ThetaDeg = ExcelRead('LOU_PN_ThetaDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOU_PN_MagDriftdB1Hz = ExcelRead('LOU_PN_MagDriftdB1Hz',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOU_PN_offset = ExcelRead('LOU_PN_offset',2,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOU_PN_f_offset_Hz = LOU_PN_offset{1};
LOU_PN_g_offset_dBc1Hz = LOU_PN_offset{2};

LOU_IMB_PhsDeg = ExcelRead('LOU_IMB_PhsDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOU_IMB_MagdB = ExcelRead('LOU_IMB_MagdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOU_SPURS = ExcelRead('LOU_SPURS',3,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOU_SPURS_foffset_spurs_Hz = LOU_SPURS{1};
LOU_SPURS_g_spurs_dBc1Hz = LOU_SPURS{2};

flagT4_LOU_PN = ExcelRead('flagT4_LOU_PN',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
flagT4_LOU_IMB = ExcelRead('flagT4_LOU_IMB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
flagT4_LOU_SPURS = ExcelRead('flagT4_LOU_SPURS',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
% flagT4_LOU_QEC = ExcelRead('flagT4_LOU_QEC',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
flagT4_LOU_QEC = 'off'

% parameter assignment
flagT4_LOU.AMtoPMPhsDriftDeg = LOU_AMtoPMPhsDriftDeg; %% 2020-03-14, g3, Add PhaseDriftDeg for AM/PM

% ********** LO Level & Frequency: **********
LOU.leveldB = LOU_leveldB;
LOU.fLO = LOU_fLO;

% ********** LO Phase Noise input: **********
LOU.PN = [];
LOU.PN.ThetaDeg = LOU_PN_ThetaDeg; % %% 2020-3-14, if PN_LOU.PN_ThetaDeg==0, no-assignment
LOU.PN.MagDriftdB1Hz = LOU_PN_MagDriftdB1Hz;

LOU.PN.f_offset_Hz = LOU_PN_f_offset_Hz; % phase noise spectrum, frequencies
LOU.PN.g_offset_dBc1Hz = LOU_PN_g_offset_dBc1Hz; % phase noise spectrum, magnitude

% ********** LO IQ Imbalance input: **********
LOU.IMB = [];
LOU.IMB.PhsDeg = LOU_IMB_PhsDeg;
LOU.IMB.MagdB = LOU_IMB_MagdB;

% ********** LO SPURS input: **********
LOU.SPURS = [];
LOU.SPURS.foffset_spurs_Hz = LOU_SPURS_foffset_spurs_Hz; % discrete spurs, freq relative to fLO
LOU.SPURS.g_spurs_dBc1Hz = LOU_SPURS_g_spurs_dBc1Hz; % discrete spurs, power relative to fLO

flagT4_LOU.PhsNoise = flagT4_LOU_PN; % LO wo PN
flagT4_LOU.IMB = flagT4_LOU_IMB; % LO wo IQ imbalance
flagT4_LOU.SPURS = flagT4_LOU_SPURS; % LO wo SPURS
flagT4_LOU.QEC = flagT4_LOU_QEC;

% ========================================

[loU_ideal,loU_realistic,tableT4_LOUInput,tableT4_LOU] = SYM_LOgenApp(LOU, fs, Nsamps, flagT4_LOU, fnum, 'half', 'LOU', []);
fnum = fnum+1

%% 1c. Generate real signal RXQECRFCSource by UpConversion + BPF
% Input 1c: ========================================
fLO = LOU_fLO 
fcuttoff_QEC = 5e6;
bwRFQEC = fix((fBBQEC+fLO+fcuttoff_QEC*[-1 1])/df)*df
flag_BPF_of_UpConversion = 'on'

% ========================================

switch flag_BPF_of_UpConversion
    case {'on'}
        [RXQECRFSource, RXQECRFSourceI, RXQECRFSourceQ] = Mixer_Up_Down_Convert_g(RXQECBBSource_real, loU_realistic(1,:), [], fs, ['RXQECRFSource'],'Up',[bwRFQEC],[fnum]);
    case {'off'}
        error('BPF is need for RXQECRFCSource')
        [RXQECRFSource, RXQECRFSourceI, RXQECRFSourceQ] = Mixer_Up_Down_Convert_g(RXQECBBSource_real, loU_realistic(1,:), [], fs, ['RXQECRFSource'],'Up',[],[fnum]);
end
fnum = fnum+1

%% 1c1. Couple the RXQECRFCSource to DownConversion Mixer, Coupler Loss and AGC Gain/NF/Unlinear
% Input 1c1, ========================================
GaindB_coupler = -30
IpwrdB_Target_coupler = []

GaindB_AGC = 0
NFdB_AGC = 4
Poip3dB_AGC = []
PDCdB_AGC = []
OP1dB_AGC = []
% ========================================

[RXQECRFSource_Coupler, evmwiGain_idC, tablewiGainOIP3_idC1] = SYM_GainIM3App(RXQECRFSource, GaindB_coupler, [], [], [], [], [], [], fs, IpwrdB_Target_coupler, bwRFQEC, [], fnum, 'Coupler', []);
fnum = fnum+1

[RXQECRFSource_Coupler_AGC, evmwiGain_idC2, tablewiGainOIP3_idC2] = SYM_GainIM3App(RXQECRFSource_Coupler, GaindB_AGC, NFdB_AGC, Poip3dB_AGC, [], PDCdB_AGC, [], OP1dB_AGC, fs, [], bwRFQEC, [], fnum, 'AGC', []);
fnum = fnum+1

%% 1d. Generate LOD for Imbalanced DownConversion Mixer
% Input 1d: ========================================
edit SystemSim_LoadCoefficient.m

LOD_AMtoPMPhsDriftDeg = ExcelRead('LOD_AMtoPMPhsDriftDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOD_leveldB = ExcelRead('LOD_leveldB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOD_fLO = ExcelRead('LOD_fLO',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOD_PN_ThetaDeg = ExcelRead('LOD_PN_ThetaDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOD_PN_MagDriftdB1Hz = ExcelRead('LOD_PN_MagDriftdB1Hz',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOD_PN_offset = ExcelRead('LOD_PN_offset',5,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOD_PN_f_offset_Hz = LOD_PN_offset{1};
LOD_PN_g_offset_dBc1Hz = LOD_PN_offset{2};

A1 = 8
LOD_IMB_PhsDeg = A1*0.5;
LOD_IMB_MagdB = A1*2;

LOD_SPURS = ExcelRead('LOD_SPURS',6,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOD_SPURS_foffset_spurs_Hz = LOD_SPURS{1}
LOD_SPURS_g_spurs_dBc1Hz = LOD_SPURS{2}

flagR2_LOD_PN = 1
flagR2_LOD_IMB = 1
flagR2_LOD_SPURS = 1;
flagR2_LOD_QEC = 'off';

% parameter assignment
flagR2_LOD.AMtoPMPhsDriftDeg = LOD_AMtoPMPhsDriftDeg; %% 2020-03-14, g3, Add PhaseDriftDeg for AM/PM

% ********** LO Level & Frequency: **********
LOD.leveldB = LOD_leveldB;
LOD.fLO = LOD_fLO;

% ********** LO Phase Noise input: **********
LOD.PN = [];
LOD.PN.ThetaDeg = LOD_PN_ThetaDeg; % %% 2020-3-14, if PN_LOU.PN_ThetaDeg==0, no-assignment
LOD.PN.MagDriftdB1Hz = LOD_PN_MagDriftdB1Hz;

LOD.PN.f_offset_Hz = LOD_PN_f_offset_Hz; % phase noise spectrum, frequencies
LOD.PN.g_offset_dBc1Hz = LOD_PN_g_offset_dBc1Hz; % phase noise spectrum, magnitude

% ********** LO IQ Imbalance input: **********
LOD.IMB = [];
LOD.IMB.PhsDeg = LOD_IMB_PhsDeg;
LOD.IMB.MagdB = LOD_IMB_MagdB;

% ********** LO SPURS input: **********
LOD.SPURS = [];
LOD.SPURS.foffset_spurs_Hz = LOD_SPURS_foffset_spurs_Hz; % discrete spurs, freq relative to fLO
LOD.SPURS.g_spurs_dBc1Hz = LOD_SPURS_g_spurs_dBc1Hz; % discrete spurs, power relative to fLO

flagR2_LOD.PhsNoise = flagR2_LOD_PN; % LO wo PN
flagR2_LOD.IMB = flagR2_LOD_IMB; % LO wo IQ imbalance
flagR2_LOD.SPURS = flagR2_LOD_SPURS; % LO wo SPURS
flagR2_LOD.QEC = flagR2_LOD_QEC;

% ========================================

[loD_ideal,loD_realistic,tableR2_LODInput,tableR2_LOD] = SYM_LOgenApp(LOD, fs, Nsamps, flagR2_LOD, fnum, 'semilogx', 'LOD', []);

%% 1e. Generate complex signal of RXQECBBCSink by Imbalanced DownConversion + w/ and w/o LPF
% Input 1e: ========================================
fcuttoff_QEC = 5e6;
bwBBQEC = fix((0+0+(fBBQEC+fcuttoff_QEC)*[-1 1])/df)*df;
flag_LPF_of_DnConversion = 'on'
% ========================================

switch flag_LPF_of_DnConversion
    case {'on'}
        [RXQECBBSink, RXQECBBCSinkI, RXQECBBCSinkQ] = Mixer_Up_Down_Convert_g(RXQECRFSource_Coupler_AGC, loD_realistic(1,:), loD_realistic(2,:), fs, ['IMB RXQECBBSink'],'Down',[bwBBQEC],[fnum]);
    case {'off'}
        [RXQECBBSink, RXQECBBCSinkI, RXQECBBCSinkQ] = Mixer_Up_Down_Convert_g(RXQECRFSource_Coupler_AGC, loD_realistic(1,:), loD_realistic(2,:), fs, ['IMB RXQECBBSink'],'Down',[],[fnum]);
end
fnum = fnum+1

%% 1f. Calculation the Imbalanced parameters from RXQECBBCSink + w/ and w/o LPF
% Input 1f: ========================================
flag_LPF_of_RXQECAnalysis = 'off'
% ========================================

switch flag_LPF_of_RXQECAnalysis
    case {'on'}
        [QECest_MagdB_RXMixer, QECest_PhsDeg_RXMixer, ~] = trx_QEC_g(RXQECBBSink, [], [], fs, [bwBBQEC], 'RXQEC','DNConv_FIR_IQDemod', 1);
    case {'off'}
        [QECest_MagdB_RXMixer, QECest_PhsDeg_RXMixer, ~] = trx_QEC_g(RXQECBBSink, [], [], fs, [], 'RXQEC','DNConv_FIR_IQDemod', 1);
end

%% 1g. Correct the RX Imbalance parameters from QECest and generate Correct loD_Corr
LOD_Corr = LOD
flagR2_LOD_Corr = flagR2_LOD

LOD_Corr.IMB.PhsDeg = LOD_Corr.IMB.PhsDeg - QECest_PhsDeg_RXMixer;
LOD_Corr.IMB.MagdB = LOD_Corr.IMB.MagdB - QECest_MagdB_RXMixer;
[loD_ideal,loD_Corr,tableR2_LODInput_Corr,tableR2_LOD_Corr] = SYM_LOgenApp(LOD_Corr, fs, Nsamps, flagR2_LOD_Corr, fnum, 'semilogx', 'LOD', []);
fnum = fnum + 1

%% 
%% 2. TX Mixer QEC *******************************************************
%% 2a. Load waveform asign for TXQECBBCSource
% Input 2a: ========================================
flag_TXQECBBSource = 'MultiCarrier'
flag_TXQECBBSource = 'SignleTone'

switch flag_TXQECBBSource
    case {'MultiCarrier'}
        
        error('QEC Not support MultiCarrier Signal!')
        
        TXQECBBSource = load('waveform_TXQECBBCSource.mat');
        TXQECBBSource = TXQECBBSource.waveformT3cell_AWGN;
        TXQECBBSource = TXQECBBSource{:}.';
        bwInband_TXQECBBSource = 9e6*[-1 1]
        bwRFQEC = fix((fLO+(20e6/2+fcuttoff_QEC)*[-1 1])/df)*df

        %% Change DownConversion fLO for MultiCarrier QEC
%         LOD_Corr. fLO = LOD.fLO - 30e6
        LOD_Corr. fLO = LOD.fLO
        [loD_ideal,loD_Corr,tableR2_LODInput_Corr,tableR2_LOD_Corr] = SYM_LOgenApp(LOD_Corr, fs, Nsamps, flagR2_LOD_Corr, fnum, 'semilogx', 'LOD', []);
        fnum = fnum + 1

    case {'SignleTone'}
        fBBQEC = 10e6
        bwRFQEC = fix((fBBQEC+fLO+fcuttoff_QEC*[-1 1])/df)*df
        TXQECBBSource = exp(1i*2*pi*fBBQEC*t);
        bwInband_TXQECBBSource = fBBQEC*[-1 1]
end
% ========================================
PLOT_FFT_dB_g(TXQECBBSource, fs, Nsamps, 'TXQECBBSource', 'df', 'full', 'pwr', fnum, [], []);
fnum = fnum+1

%% 2b. Generate LOU for Imbalanced UpConversion Mixer
% Input 2b: ========================================

edit SystemSim_LoadCoefficient.m

LOU_AMtoPMPhsDriftDeg = ExcelRead('LOU_AMtoPMPhsDriftDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOU_leveldB = ExcelRead('LOU_leveldB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOU_fLO = ExcelRead('LOU_fLO',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOU_PN_ThetaDeg = ExcelRead('LOU_PN_ThetaDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOU_PN_MagDriftdB1Hz = ExcelRead('LOU_PN_MagDriftdB1Hz',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOU_PN_offset = ExcelRead('LOU_PN_offset',2,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOU_PN_f_offset_Hz = LOU_PN_offset{1};
LOU_PN_g_offset_dBc1Hz = LOU_PN_offset{2};

LOU_IMB_PhsDeg = ExcelRead('LOU_IMB_PhsDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOU_IMB_MagdB = ExcelRead('LOU_IMB_MagdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

LOU_SPURS = ExcelRead('LOU_SPURS',3,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
LOU_SPURS_foffset_spurs_Hz = LOU_SPURS{1};
LOU_SPURS_g_spurs_dBc1Hz = LOU_SPURS{2};

flagT4_LOU_PN = ExcelRead('flagT4_LOU_PN',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
flagT4_LOU_IMB = ExcelRead('flagT4_LOU_IMB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
flagT4_LOU_SPURS = ExcelRead('flagT4_LOU_SPURS',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
flagT4_LOU_QEC = ExcelRead('flagT4_LOU_QEC',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);

% parameter assignment
flagT4_LOU.AMtoPMPhsDriftDeg = LOU_AMtoPMPhsDriftDeg; %% 2020-03-14, g3, Add PhaseDriftDeg for AM/PM

% ********** LO Level & Frequency: **********
LOU.leveldB = LOU_leveldB;
LOU.fLO = LOU_fLO;

% ********** LO Phase Noise input: **********
LOU.PN = [];
LOU.PN.ThetaDeg = LOU_PN_ThetaDeg; % %% 2020-3-14, if PN_LOU.PN_ThetaDeg==0, no-assignment
LOU.PN.MagDriftdB1Hz = LOU_PN_MagDriftdB1Hz;

LOU.PN.f_offset_Hz = LOU_PN_f_offset_Hz; % phase noise spectrum, frequencies
LOU.PN.g_offset_dBc1Hz = LOU_PN_g_offset_dBc1Hz; % phase noise spectrum, magnitude

% ********** LO IQ Imbalance input: **********
LOU.IMB = [];
LOU.IMB.PhsDeg = LOU_IMB_PhsDeg;
LOU.IMB.MagdB = LOU_IMB_MagdB;

% ********** LO SPURS input: **********
LOU.SPURS = [];
LOU.SPURS.foffset_spurs_Hz = LOU_SPURS_foffset_spurs_Hz; % discrete spurs, freq relative to fLO
LOU.SPURS.g_spurs_dBc1Hz = LOU_SPURS_g_spurs_dBc1Hz; % discrete spurs, power relative to fLO


flagT4_LOU.PhsNoise = flagT4_LOU_PN; % LO wo PN
flagT4_LOU.IMB = 1; % LO wo IQ imbalance
flagT4_LOU.SPURS = flagT4_LOU_SPURS; % LO wo SPURS
flagT4_LOU.QEC = 'off';

if flagT4_LOU_PN==0
    [loU_ideal,loU_realistic,tableT4_LOUInput,tableT4_LOU] = SYM_LOgenApp(LOU, fs, Nsamps, flagT4_LOU, fnum, 'half', 'LOU', []);
else
    [loU_ideal,loU_realistic,tableT4_LOUInput,tableT4_LOU] = SYM_LOgenApp(LOU, fs, Nsamps, flagT4_LOU, fnum, 'semilogx', 'LOU', []);
end

flagT4_LOU_Unlinearity = ExcelRead('flagT4_LOU_UNL',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
flagT4_LOU_Unlinearity = 'off'

if ~strcmp(flagT4_LOU_Unlinearity,'off')
    % input
    LOU_GaindB = ExcelRead('LOU_GaindB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_NFdB = ExcelRead('LOU_NFdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_Poip3dB = ExcelRead('LOU_Poip3dB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_PHD2dB = ExcelRead('LOU_PHD2dB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_PDCdB = ExcelRead('LOU_PDCdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_AMtoPMDegDrift = ExcelRead('LOU_AMtoPMDegDrift',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_OP1dB = ExcelRead('LOU_OP1dB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_IpwrdB_Target = ExcelRead('LOU_IpwrdB_Target',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_bwACLROffset = ExcelRead('LOU_bwACLROffset',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    fnum = fnum+1;
    [loU_realistic_UNL, ~,  tableT4_LOU_UNL] = SYM_GainIM3App(loU_realistic, LOU_GaindB, LOU_NFdB, LOU_Poip3dB, LOU_PHD2dB, LOU_PDCdB, LOU_AMtoPMDegDrift, LOU_OP1dB, fs, LOU_IpwrdB_Target, [], LOU_bwACLROffset, fnum, {flagT4_LOU_Unlinearity}, fnum_dir);
    % export
    loU_output = loU_realistic_UNL;
else
    % export
    loU_output = loU_realistic;
end

%% 2c. Generate TXQECRFSource by Imbalanced UpConversion + w/ BPF
flag_BPF_of_UpConversion_TXQEC = 'on'

switch flag_BPF_of_UpConversion_TXQEC
    case {'on'}
        [TXQECRFSource, TXQECRFSourceI, TXQECRFSourceQ] = Mixer_Up_Down_Convert_g(TXQECBBSource, loU_output(1,:), loU_output(2,:), fs, ['TXQECRFSource'], 'Up',[bwRFQEC],fnum);
    case {'off'}
        error('BPF for Image Rejection is need for TXQECRFSource, or the QEC phase will be failure!')
end
fnum = fnum+1

%% 2d. Generate complex signal of TXQECBBSink by DownConversion + w/ and w/o LPF
% Input 2d: ========================================
bwBBQEC = fix((mean(bwRFQEC) - LOD_Corr.fLO + 20e6/2+fcuttoff_QEC)*[-1 1]/df)*df
flag_LPF_of_DnConversion_TXQEC = 'on'
% ========================================

switch flag_LPF_of_DnConversion_TXQEC
    case {'on'}
%         [TXQECBBSink, TXQECBBSinkI, TXQECBBSinkQ] = Mixer_Up_Down_Convert_g(TXQECRFSource, loD_Corr(1,:), loD_Corr(2,:), fs, ['IMB TXQECBBSink'],'Down',[bwBBQEC],[fnum]);
        [TXQECBBSink, TXQECBBSinkI, TXQECBBSinkQ] = Mixer_Up_Down_Convert_g([TXQECRFSourceI;TXQECRFSourceQ], loD_Corr(1,:), loD_Corr(2,:), fs, ['IMB TXQECBBSink'],'Down',[bwBBQEC],[fnum]);
    case {'off'}
        [TXQECBBSink, TXQECBBSinkI, TXQECBBSinkQ] = Mixer_Up_Down_Convert_g([TXQECRFSourceI;TXQECRFSourceQ], loD_Corr(1,:), loD_Corr(2,:), fs, ['IMB TXQECBBSink'],'Down',[],[fnum]);
end
fnum = fnum+1

%% 2e. Calculation the Imbalanced parameters from TXQECBBSink + w/ LPF
% Input 2e: ========================================
if strcmp(flag_LPF_of_DnConversion_TXQEC,'off')
    flag_LPF_of_TXQECAnalysis = 'on'
elseif strcmp(flag_LPF_of_DnConversion_TXQEC,'on')
    flag_LPF_of_TXQECAnalysis = 'off'
else
    error('Add LPF before RXQEC!')
end
% ========================================

switch flag_LPF_of_TXQECAnalysis
    case {'on'}
%         [QECest_MagdB_TXMixer, QECest_PhsDeg_TXMixer, ~] = trx_QEC_g(TXQECBBSink, [], [], fs, bwBBQEC, 'TXQEC','DNConv_FIR_IQDemod', 1);
        [QECest_MagdB_TXMixer, QECest_PhsDeg_TXMixer, ~] = trx_QEC_g([TXQECBBSinkI;TXQECBBSinkQ], [], [], fs, bwBBQEC, 'RXQEC','DNConv_FIR_IQDemod', 1);
%         [QECest_MagdB_TXMixer, QECest_PhsDeg_TXMixer, ~] = trx_QEC_g([TXQECBBSinkI;TXQECBBSinkQ], [], [], fs, bwBBQEC, 'TXQEC','DNConv_FIR_IQDemod', 1);
    case {'off'}
        warning('LPF is need for TXQECBBSink')
        [QECest_MagdB_TXMixer, QECest_PhsDeg_TXMixer, ~] = trx_QEC_g([TXQECBBSinkI;TXQECBBSinkQ], [], [], fs, [], 'RXQEC','DNConv_FIR_IQDemod', 1);
end
