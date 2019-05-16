subjects = [];
valid_subjects=2:9;
%% pilot
% subjects(1).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/pilot/matspilot1_1200hz_20150330_01.ds';
% subjects(1).trialfun = 'mytrialfun_pilot';
projectdir = '/project/3011085.02/phasecode/';
datadir = [projectdir, 'DATA/'];
%% experiment
% 1
subjects(1).sessions = [1,2,3];
subjects(1).channels = {'MEG','-MLF14', '-MLF32', '-MRO51', '-MLT31', '-MLT41','-MRT31', '-MRT41', '-MRT51'}; %, '-MRT11', '-MRT12'};
subjects(1).session(1).dataset = [datadir, '3016045.07_matves_001_001/20150609matves01_1200hz_20150609_01.ds'];
subjects(1).session(1).icacomp = [7 29 36];
% subjects(1).session(1).icacompLP = []; %only a jump, no ecg
subjects(1).session(2).dataset = [datadir, '3016045.07_matves_001_002/matves01b_1200hz_20150610_01.ds'];
subjects(1).session(2).icacomp = [7 33 50];
% subjects(1).session(2).icacompLP = [];%no ecg
subjects(1).session(3).dataset = [datadir, '3016045.07_matves_001_003/301604507matves013_1200hz_20150611_01.ds'];
subjects(1).session(3).icacomp = [15 35];
% subjects(1).session(3).icacompLP = [];
subjects(1).trialfun = 'mytrialfun';
subjects(1).class{1} = 328;
subjects(1).class{2}.hemifield(1) = 241; %ntrials
subjects(1).class{2}.hemifield(2) = 261;
subjects(1).class{3} = 517;

% 2

subjects(2).sessions = [1,2,3];
subjects(2).channels = {'MEG','-MLT31', '-MLT41', '-MLF14', '-MRT41', '-MRT42', '-MRT51', '-MRT52'};
subjects(2).session(1).dataset = [datadir, '3016045.07_matves_002_001/301604507matves021_1200hz_20150622_01.ds'];
subjects(2).session(1).icacomp = [20];
% subjects(2).session(1).icacompLP = [];
subjects(2).session(2).dataset = [datadir, '3016045.07_matves_002_002/301604507matves022_1200hz_20150623_01.ds'];
subjects(2).session(2).icacomp = [26];
% subjects(2).session(2).icacompLP = [];
subjects(2).session(3).dataset = [datadir, '3016045.07_matves_002_003/301604507matves023_1200hz_20150624_01.ds'];
subjects(2).session(3).icacomp = [15 37];
% subjects(2).session(3).icacompLP = [];
subjects(2).trialfun = 'mytrialfun';
subjects(2).class{1} = 318;
subjects(2).class{2}.hemifield(1) = 219;
subjects(2).class{2}.hemifield(2) = 285;
subjects(2).class{3} = 451;

% 3
subjects(3).sessions = [1,12,2,3];
subjects(3).channels = {'MEG','-MLF14', '-MLF32','-MLT31', '-MRT31', '-MRT41', '-MRT56'};
subjects(3).session(1).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject03/session1/301604507matves031_1200hz_20150624_01.ds';
subjects(3).session(1).icacomp = [10 13];
% subjects(3).session(1).icacompLP = [];
subjects(3).session(12).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject03/session1/301604507matves031b_1200hz_20150624_01.ds';%session 1b
subjects(3).session(12).icacomp = [5 11 21];
% subjects(3).session(12).icacompLP = [];
subjects(3).session(2).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject03/session2/301604507matves032_1200hz_20150625_01.ds';
subjects(3).session(2).icacomp = [10 23];
% subjects(3).session(2).icacompLP = [];
subjects(3).session(3).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject03/session3/301604507matves033_1200hz_20150626_01.ds';
subjects(3).session(3).icacomp = [8 30];
% subjects(3).session(3).icacompLP = [];
subjects(3).trialfun = 'mytrialfun';
subjects(3).class{1}=249;
subjects(3).class{2}.hemifield(1) = 203;
subjects(3).class{2}.hemifield(2) = 223;
subjects(3).class{3} = 435;

% 4
subjects(4).sessions = [1,2,22,3];
subjects(4).channels = {'MEG', '-MLT31', '-MLT41', '-MLF32'};
subjects(4).session(1).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject04/session1/301604507matves041_1200hz_20150709_01.ds';
subjects(4).session(1).icacomp = [6 12 32];
% subjects(4).session(1).icacompLP = [];
subjects(4).session(2).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject04/session2/301604507matves042_1200hz_20150710_01.ds';
subjects(4).session(2).icacomp = [6 13 18];
% subjects(4).session(2).icacompLP = [];
subjects(4).session(22).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject04/session2/301604507matves042b_1200hz_20150710_01.ds';
subjects(4).session(22).icacomp = [8 13];
% subjects(4).session(22).icacompLP = [];
subjects(4).session(3).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject04/session3/301604507matves043_1200hz_20150713_01.ds';
subjects(4).session(3).icacomp = [5 42];
% subjects(4).session(3).icacompLP = [];
subjects(4).trialfun = 'mytrialfun';
subjects(4).class{1}=186;
subjects(4).class{2}.hemifield(1) = 132;
subjects(4).class{2}.hemifield(2) = 150;
subjects(4).class{3} = 275;

% 5
subjects(5).sessions = [1,2,3];
subjects(5).channels = {'MEG', '-MRT31'};
subjects(5).session(1).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject05/session1/301604507matves051_1200hz_20150713_01.ds';
subjects(5).session(1).icacomp = [12 30];
% subjects(5).session(1).icacompLP = [];
subjects(5).session(2).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject05/session2/301604507matves052_1200hz_20150714_01.ds';
subjects(5).session(2).icacomp = [9 21];
% subjects(5).session(2).icacompLP = [];
subjects(5).session(3).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject05/session3/301604507matves053_1200hz_20150715_01.ds';
subjects(5).session(3).icacomp = [12 28];
% subjects(5).session(3).icacompLP = [];
subjects(5).trialfun = 'mytrialfun';
subjects(5).class{1} = 395;
subjects(5).class{2}.hemifield(1) = 329;
subjects(5).class{2}.hemifield(2) = 320;
subjects(5).class{3} = 660;

% 6
subjects(6).sessions = [1,2,3];
subjects(6).channels = {'MEG','-MLF32'};
subjects(6).session(1).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject06/session1/301604507matves061_1200hz_20150720_01.ds';
subjects(6).session(1).icacomp = [12 19];
% subjects(6).session(1).icacompLP = [];
subjects(6).session(2).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject06/session2/301604507matves062_1200hz_20150722_01.ds';
subjects(6).session(2).icacomp = [10 13 17 ];
% subjects(6).session(2).icacompLP = [];
subjects(6).session(3).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject06/session3/301604507matves063_1200hz_20150723_01.ds';
subjects(6).session(3).icacomp = [3 9 33 ];
% subjects(6).session(3).icacompLP = [];
subjects(6).trialfun = 'mytrialfun';
subjects(6).class{1} = 385;
subjects(6).class{2}.hemifield(1) = 299;
subjects(6).class{2}.hemifield(2) = 326;
subjects(6).class{3} = 600;

% 7
subjects(7).sessions = [1,2,3];
subjects(7).channels = {'MEG', '-MRT31', '-MLF32', '-MRT41'};
subjects(7).session(1).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject07/session1/301604507matves071_1200hz_20150720_01.ds';
subjects(7).session(1).icacomp = [11];
% subjects(7).session(1).icacompLP = [];
subjects(7).session(2).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject07/session2/301604507matves072_1200hz_20150722_01.ds';
subjects(7).session(2).icacomp = [7 11];
% subjects(7).session(2).icacompLP = [];
subjects(7).session(3).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject07/session3/301604507matves073_1200hz_20150723_01.ds';
subjects(7).session(3).icacomp = [19];
% subjects(7).session(3).icacompLP = [];
subjects(7).trialfun = 'mytrialfun';
subjects(7).class{1} = 333;
subjects(7).class{2}.hemifield(1) = 232;
subjects(7).class{2}.hemifield(2) = 279;
subjects(7).class{3} = 501;

% 8
subjects(8).sessions = [1,2,3];
subjects(8).channels = {'MEG', '-MRT31'};
subjects(8).session(1).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject08/session1/301604507matves081_1200hz_20150727_01.ds';
subjects(8).session(1).icacomp = [15 37 ];
% subjects(8).session(1).icacompLP = [];
subjects(8).session(2).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject08/session2/301604507matves082_1200hz_20150728_01.ds';
subjects(8).session(2).icacomp = [23];
% subjects(8).session(2).icacompLP = [];
subjects(8).session(3).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject08/session3/301604507matves083_1200hz_20150729_01.ds';
subjects(8).session(3).icacomp = [2 22];
% subjects(8).session(3).icacompLP = [];
subjects(8).trialfun = 'mytrialfun';
subjects(8).class{1}=324;
subjects(8).class{2}.hemifield(1) = 294;
subjects(8).class{2}.hemifield(2) = 260;
subjects(8).class{3} = 538;

% 9
subjects(9).sessions = [1,2,3,32];
subjects(9).channels = {'MEG', '-MLF12'};
subjects(9).session(1).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject09/session1/301604507matves091_1200hz_20150727_01.ds';
subjects(9).session(1).icacomp = [15 19];
% subjects(9).session(1).icacompLP = [];
subjects(9).session(2).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject09/session2/301604507matves092_1200hz_20150728_01.ds';
subjects(9).session(2).icacomp = [8 17];
% subjects(9).session(2).icacompLP = [];
subjects(9).session(3).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject09/session3/301604507matves093_1200hz_20150729_01.ds';
subjects(9).session(3).icacomp = [10 12];
% subjects(9).session(3).icacompLP = [];
subjects(9).session(32).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject09/session3/301604507matves093b_1200hz_20150729_01.ds';
subjects(9).session(32).icacomp = [9 19];
% subjects(9).session(32).icacompLP = [];
subjects(9).trialfun = 'mytrialfun';
subjects(9).class{1} = 350;
subjects(9).class{2}.hemifield(1) = 245;
subjects(9).class{2}.hemifield(2) = 302;
subjects(9).class{3} = 520;

% 10
subjects(10).sessions = [1,2,22,3];
subjects(10).channels = {'MEG'};
subjects(10).session(1).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject10/session1/301604507matves101_1200hz_20150810_01.ds';
subjects(10).session(1).icacomp = [5 7 14];
% subjects(10).session(1).icacompLP = [];
subjects(10).session(2).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject10/session2/301604507matves102_1200hz_20150811_01.ds';
subjects(10).session(2).icacomp = [3 5];
% subjects(10).session(2).icacompLP = [];
subjects(10).session(22).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject10/session2/301604507matves102c_1200hz_20150811_01.ds';
subjects(10).session(22).icacomp = [5 15];
% subjects(10).session(22).icacompLP = [];
subjects(10).session(3).dataset = '/home/electromag/matves/Data/phasecode/meg_raw/subject10/session3/301604507matves103_1200hz_20150812_01.ds';
subjects(10).session(3).icacomp = [15 22];
% subjects(10).session(3).icacompLP = [];
subjects(10).trialfun = 'mytrialfun';
subjects(10).class{1} = 227;
subjects(10).class{2}.hemifield(1) = 194;
subjects(10).class{2}.hemifield(2) = 178;
subjects(10).class{3} = 359;










%% Number of rejected trials
% 01: 163, 199
% 