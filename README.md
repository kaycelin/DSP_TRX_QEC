# DSP_QEC
Background: 
- For Zero_IF architecture the mixer imbalance will introduce the image and damage the SNR (EVM)
- To calculation the IQ imbalance and compensation by tuning the hardware IQ-mixer OR compensation by DSP

Experiment:
- Design LOI and LOQ with magnitude and phase imbalance
- Simulate the hardware behavior and calcuation the IQ imbalance by QEC Receiver

Hardware architecture:
- Upconversion mixer and Downconversion mixer
- QEC source (CW or Multicarrier) and QEC source feedback path
- QEC receiver with ADC for DSP to calculate the IQ imbalance
- Hardward IQ-mixer tunning control path 
![image](https://user-images.githubusercontent.com/87049112/142336782-b250713c-13d8-45c3-8e4a-34bda1699515.png)

Process: 
# RX QEC: 
1a. Generate real signal of RXQECBBCSource_real, 10MHz          
- Input 1a:
  - fs = 3.9322e+09
  - Nsamps = 2457600
  - df = fs/Nsamps
  - fBBQEC = 10e6
  - fnum = 111801

1b. Generate LOU for Real Signal UpConversion
- Input 1b:
  - LOU_AMtoPMPhsDriftDeg = 0
  - LOU_leveldB = 20
  - LOU_fLO = 300e6

  - LOU_PN_ThetaDeg = 1
  - LOU_PN_MagDriftdB1Hz = 1

  - LOU_PN_f_offset_Hz = 1e9*[ 0    0.0001    0.0010    0.0050    0.0100    0.0500    0.1000    1.0000 ]
  - LOU_PN_g_offset_dBc1Hz = 1e9*[ 0   -80  -100  -120  -130  -130  -130  -130 ]

  - LOU_IMB_PhsDeg = -12
  - LOU_IMB_MagdB = 8

  - LOU_SPURS_foffset_spurs_Hz = 1e3[ 10         100        1000        1500 ]
  - LOU_SPURS_g_spurs_dBc1Hz = [ -60   -60   -60   -60 ]

  - flagT4_LOU_PN = 0
  - flagT4_LOU_IMB = 1
  - flagT4_LOU_SPURS = 0
  - flagT4_LOU_QEC = 'off'

1c. Generate real signal RXQECRFCSource by UpConversion + BPF
- Input 1c: 
  - fLO = LOU_fLO 
  - fcuttoff_QEC = 5e6;
  - bwRFQEC = fix((fBBQEC+fLO+fcuttoff_QEC*[-1 1])/df)*df
  - flag_BPF_of_UpConversion = 'on'
![image](https://user-images.githubusercontent.com/87049112/142348291-d0a5f5e1-3890-4664-8664-11f873302848.png)

1c1. Couple the RXQECRFCSource to DownConversion Mixer, Coupler Loss and AGC Gain/NF/Unlinear
- Input 1c1:
  - GaindB_coupler = -30
  - IpwrdB_Target_coupler = []

  - GaindB_AGC = 0
  - NFdB_AGC = 4
  - Poip3dB_AGC = []
  - PDCdB_AGC = []
  - OP1dB_AGC = []
  - 
1d. Generate LOD for Imbalanced DownConversion Mixer
- Input 1d:

  - LOD_AMtoPMPhsDriftDeg = 0
  - LOD_leveldB = 20
  - LOD_fLO = 300e6

  - LOD_PN_ThetaDeg = 1
  - LOD_PN_MagDriftdB1Hz = 1

   -LOD_PN_offset = ExcelRead('LOD_PN_offset',5,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
   -LOD_PN_f_offset_Hz = 1e3*[0         100         200         400         600         800        1200        1800        6000        10000]
   -LOD_PN_g_offset_dBc1Hz = [0   -80   -95  -100  -109  -112  -115  -120  -130  -133];

  - LOD_IMB_PhsDeg = 4
  - LOD_IMB_MagdB = 16

  - LOD_SPURS_foffset_spurs_Hz = 1e3*[50         100         400       10000]
  - LOD_SPURS_g_spurs_dBc1Hz = [-30   -40   -50   -70]

  - flagR2_LOD_PN = 1
  - flagR2_LOD_IMB = 1
  - flagR2_LOD_SPURS = 1;
  - flagR2_LOD_QEC = 'off'

1e. Generate complex signal of RXQECBBCSink by Imbalanced DownConversion + w/ and w/o LPF
- Input 1e: 
  - fcuttoff_QEC = 5e6;
  - bwBBQEC = fix((0+0+(fBBQEC+fcuttoff_QEC)*[-1 1])/df)*df;
  - flag_LPF_of_DnConversion = 'on'
 ![image](https://user-images.githubusercontent.com/87049112/142348647-a240de5f-aa4c-4272-a0e4-2cd252624b46.png)

- Output Result: calculation error<0.2
  - QECest_MagdB_RXMixer = 15.9858 (compare to LOD_IMB_MagdB=16 calculation error<0.2)
  - QECest_PhsDeg_RXMixer = 3.9516 (compare to LOD_IMB_PhsDeg=4 calculation error<0.2)
 ![image](https://user-images.githubusercontent.com/87049112/142348703-4ffc36b8-bc2c-4a10-9bbb-0a1a5bc079a3.png)

1g. Correct the RX Imbalance parameters from QECest and generate Correct loD_Corr

# TX QEC
2a. Load waveform asign for TXQECBBCSource
- Input 2a:
  - flag_TXQECBBSource = 'SignleTone'; (QEC Not support MultiCarrier Signal?)

2b. Generate LOU for Imbalanced UpConversion Mixer
- Input 2b: (the same Input 1b)

  - LOU_AMtoPMPhsDriftDeg = 0
  - LOU_leveldB = 20
  - LOU_fLO = 300e6

  - LOU_PN_ThetaDeg = 1
  - LOU_PN_MagDriftdB1Hz = 1

  - LOU_PN_f_offset_Hz = 1e9*[ 0    0.0001    0.0010    0.0050    0.0100    0.0500    0.1000    1.0000 ]
  - LOU_PN_g_offset_dBc1Hz = 1e9*[ 0   -80  -100  -120  -130  -130  -130  -130 ]

  - LOU_IMB_PhsDeg = -12
  - LOU_IMB_MagdB = 8

  - LOU_SPURS_foffset_spurs_Hz = 1e3[ 10         100        1000        1500 ]
  - LOU_SPURS_g_spurs_dBc1Hz = [ -60   -60   -60   -60 ]

  - flagT4_LOU_PN = 0
  - flagT4_LOU_IMB = 1
  - flagT4_LOU_SPURS = 0
  - flagT4_LOU_QEC = 'off'

2c. Generate TXQECRFSource by Imbalanced UpConversion + w/ BPF
  - flag_BPF_of_UpConversion_TXQEC = 'on'; (BPF for Image Rejection is need for TXQECRFSource, or the QEC phase will be failure!)
![image](https://user-images.githubusercontent.com/87049112/142348994-e3d85784-2210-4947-b490-8a65122bb1a6.png)

2d. Generate complex signal of TXQECBBSink by DownConversion + w/ and w/o LPF
- Input 2d: ========================================
  - bwBBQEC = fix((mean(bwRFQEC) - LOD_Corr.fLO + 20e6/2+fcuttoff_QEC)*[-1 1]/df)*df
  - flag_LPF_of_DnConversion_TXQEC = 'on'

2e. Calculation the Imbalanced parameters from TXQECBBSink + w/ LPF
- Output Result: calculation error<0.1
  - QECest_MagdB_TXMixer = 8.0016 (compare to LOU_IMB_MagdB=8 calculation error<0.1)
  - QECest_PhsDeg_TXMixer = -12.0501 (compare to LOU_IMB_PhsDeg=-12 calculation error<0.1)
![image](https://user-images.githubusercontent.com/87049112/142349162-ce173054-6f00-443a-825c-dc39ddb13893.png)





