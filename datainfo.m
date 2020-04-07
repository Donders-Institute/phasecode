subjects = [];
valid_subjects=1:10;
%% directories
projectdir = '/project/3011085.06/';
datadir = [projectdir, 'data/'];
figures_dir = [projectdir, 'manuscript/figures/'];
results_dir = [projectdir, 'results/'];

%% subject info
% subjects(1).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/pilot/matspilot1_1200hz_20150330_01.ds';
% subjects(1).trialfun = 'mytrialfun_pilot';

% 1
subjects(1).sessions = [1,2,3];
subjects(1).validsessions = [1,2,3]; % of the preprocessed sessions.
subjects(1).channels = {'MEG','-MLF14', '-MLF32', '-MRO51', '-MLT31', '-MLT41','-MRT31', '-MRT41', '-MRT51'}; %, '-MRT11', '-MRT12'};
subjects(1).session(1).dataset = [datadir, 'sub01/meg01/20150609matves01_1200hz_20150609_01.ds'];
subjects(1).session(2).dataset = [datadir, 'sub01/meg02/matves01b_1200hz_20150610_01.ds'];
subjects(1).session(3).dataset = [datadir, 'sub01/meg03/301604507matves013_1200hz_20150611_01.ds'];
subjects(1).mri = [datadir, '/sub01/mri/201506012_ELESTJ_S01.MR.DCCN_SEQUENCES_STANDARD_SEQUENCES.0002.0001.2015.06.12.10.29.33.349566.91157570.IMA'];
subjects(1).trialfun = 'mytrialfun';

% 2
subjects(2).sessions = [1,2,3];
subjects(2).validsessions = [1,2,3];
subjects(2).channels = {'MEG','-MLT31', '-MLT41', '-MLF14', '-MRT41', '-MRT42', '-MRT51', '-MRT52'};
subjects(2).session(1).dataset = [datadir, 'sub02/meg01/301604507matves021_1200hz_20150622_01.ds'];
subjects(2).session(2).dataset = [datadir, 'sub02/meg02/301604507matves022_1200hz_20150623_01.ds'];
subjects(2).session(3).dataset = [datadir, 'sub02/meg03/301604507matves023_1200hz_20150624_01.ds'];
subjects(2).trialfun = 'mytrialfun';
subjects(2).mri = [datadir, '/sub02/mri/TOMMAR_EELSPA_SACCADE_77704.MR.TOMMAR_SKYRA.0007.0001.2014.09.11.09.38.14.448977.370847816.IMA'];

% 3
subjects(3).sessions = [1,12,2,3];
subjects(3).validsessions = [1,2,3];
subjects(3).channels = {'MEG','-MLF14', '-MLF32','-MLT31', '-MRT31', '-MRT41', '-MRT56'};
subjects(3).session(1).dataset = [datadir, 'sub03/meg01/301604507matves031_1200hz_20150624_01.ds'];
subjects(3).session(12).dataset = [datadir, 'sub03/meg01/301604507matves031b_1200hz_20150624_01.ds'];%session 1b
subjects(3).session(2).dataset = [datadir, 'sub03/meg02/301604507matves032_1200hz_20150625_01.ds'];
subjects(3).session(3).dataset = [datadir, 'sub03/meg03/301604507matves033_1200hz_20150626_01.ds'];
subjects(3).trialfun = 'mytrialfun';
subjects(3).mri = [datadir, '/sub03/mri/ELESTJ_20141107_S32.MR.DCCN_SEQUENCES_STANDARD_SEQUENCES.0002.0001.2014.11.07.16.51.27.247752.49776744.IMA'];

% 4
subjects(4).sessions = [1,2,22,3];
subjects(4).validsessions = [2,3]; % the only one with one invalid session (1 session with chance level behavior)
subjects(4).channels = {'MEG', '-MLT31', '-MLT41', '-MLF32'};
subjects(4).session(1).dataset = [datadir, 'sub04/meg01/301604507matves041_1200hz_20150709_01.ds'];
subjects(4).session(2).dataset = [datadir, 'sub04/meg02/301604507matves042_1200hz_20150710_01.ds'];
subjects(4).session(22).dataset = [datadir, 'sub04/meg02/301604507matves042b_1200hz_20150710_01.ds'];
subjects(4).session(3).dataset = [datadir, 'sub04/meg03/301604507matves043_1200hz_20150713_01.ds'];
subjects(4).mri = [datadir, 'sub04/mri/EVGBED_20150630_P2S1.MR.EVGBED_PRISMA.0013.0001.2015.06.30.14.58.35.53572.792031346.IMA'];
subjects(4).trialfun = 'mytrialfun';

% 5
subjects(5).sessions = [1,2,3];
subjects(5).validsessions = [1,2,3];
subjects(5).session(1).dataset = [datadir, 'sub05/meg01/301604507matves051_1200hz_20150713_01.ds'];
subjects(5).session(2).dataset = [datadir, 'sub05/meg02/301604507matves052_1200hz_20150714_01.ds'];
subjects(5).session(3).dataset = [datadir, 'sub05/meg03/301604507matves053_1200hz_20150715_01.ds'];
subjects(5).trialfun = 'mytrialfun';
subjects(5).mri = [datadir, 'sub05/mri/P027_t1.nii.gz'];

% 6
subjects(6).sessions = [1,2,3];
subjects(6).validsessions = [1,2,3];
subjects(6).channels = {'MEG','-MLF32'};
subjects(6).session(1).dataset = [datadir, 'sub06/meg01/301604507matves061_1200hz_20150720_01.ds'];
subjects(6).session(2).dataset = [datadir, 'sub06/meg02/301604507matves062_1200hz_20150722_01.ds'];
subjects(6).session(3).dataset = [datadir, 'sub06/meg03/301604507matves063_1200hz_20150723_01.ds'];
subjects(6).trialfun = 'mytrialfun';
subjects(6).mri = [datadir, 'sub06/mri/P018_t1.nii.gz'];

% 7
subjects(7).sessions = [1,2,3];
subjects(7).validsessions = [1,2,3];
subjects(7).channels = {'MEG', '-MRT31', '-MLF32', '-MRT41'};
subjects(7).session(1).dataset = [datadir, 'sub07/meg01/301604507matves071_1200hz_20150720_01.ds'];
subjects(7).session(2).dataset = [datadir, 'sub07/meg02/301604507matves072_1200hz_20150722_01.ds'];
subjects(7).session(3).dataset = [datadir, 'sub07/meg03/301604507matves073_1200hz_20150723_01.ds'];
subjects(7).trialfun = 'mytrialfun';
subjects(7).mri = [datadir, 'sub07/mri/PAUGAA_20150723_106312.MR.DCCN_SEQUENCES_STANDARD_SEQUENCES.0002.0001.2015.07.23.16.37.32.976756.102123638.IMA'];

% 8
subjects(8).sessions = [1,2,3];
subjects(8).validsessions = [1,2,3];
subjects(8).channels = {'MEG', '-MRT31'};
subjects(8).session(1).dataset = [datadir, 'sub08/meg01/301604507matves081_1200hz_20150727_01.ds'];
subjects(8).session(2).dataset = [datadir, 'sub08/meg02/301604507matves082_1200hz_20150728_01.ds'];
subjects(8).session(3).dataset = [datadir, 'sub08/meg03/301604507matves083_1200hz_20150729_01.ds'];
subjects(8).trialfun = 'mytrialfun';
subjects(8).mri = [datadir, 'sub08/mri/LIEVLIE_20150616_S92.MR.LIEVLIE_FEATURE_EXP.0002.0001.2015.06.16.16.13.03.531642.240004270.IMA'];

% 9
subjects(9).sessions = [1,2,3,32];
subjects(9).validsessions = [1,2,3];
subjects(9).channels = {'MEG', '-MLF12'};
subjects(9).session(1).dataset = [datadir, 'sub09/meg01/301604507matves091_1200hz_20150727_01.ds'];
subjects(9).session(2).dataset = [datadir, 'sub09/meg02/301604507matves092_1200hz_20150728_01.ds'];
subjects(9).session(3).dataset = [datadir, 'sub09/meg03/301604507matves093_1200hz_20150729_01.ds'];
subjects(9).session(32).dataset = [datadir, 'sub09/meg03/301604507matves093b_1200hz_20150729_01.ds'];
subjects(9).trialfun = 'mytrialfun';
subjects(9).mri = [datadir, 'sub09/mri/EGBHAR_20150326_PP19.MR.EGBHAR_AVANTO.0002.0001.2015.03.26.13.01.09.48022.71055444.IMA'];

% 10
subjects(10).sessions = [1,2,22,3];
subjects(10).validsessions = [1,2,3];
subjects(10).channels = {'MEG'};
subjects(10).session(1).dataset = [datadir, 'sub10/meg01/301604507matves101_1200hz_20150810_01.ds'];
subjects(10).session(2).dataset = [datadir, 'sub10/meg02/301604507matves102_1200hz_20150811_01.ds'];
subjects(10).session(22).dataset = [datadir, 'sub10/meg02/301604507matves102c_1200hz_20150811_01.ds'];
subjects(10).session(3).dataset = [datadir, 'sub10/meg03/301604507matves103_1200hz_20150812_01.ds'];
subjects(10).trialfun = 'mytrialfun';
subjects(10).mri = [datadir, 'sub10/mri/PAUGAA_20141110_S04.MR.DCCN_SEQUENCES_STANDARD_SEQUENCES.0002.0001.2014.11.10.16.38.38.18220.50186354.IMA'];

%% variables
useparc = {'_7_B05_08','_7_B05_04','_7_B05_09','_7_B05_10','_7_B05_11',...
  '_7_B05_12','_7_B05_13', '_19_B05_12', '_19_B05_09', '_8_B05_06', ...
  '_7_B05_01','_7_B05_02','_7_B05_05','_7_B05_09', '_8_B05_06', ...
  '_5_B05_02', '_5_B05_01', '_2_B05_08','_8_B05_02','_8_B05_03'};

