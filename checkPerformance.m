%% 1st performance check: Test conditions and keyPress
% cd('/home/electromag/matves/Data/phasecode/behavioral');
cd('M:/Data/phasecode/behavioral/subject10');
load('matves10a20150810');
% load('test.mat');
numberOfBlocks = 10;
trialsPerBlock = 80;
lastTrial = numberOfBlocks * trialsPerBlock;
presented = conditions(1:lastTrial,2);
response = keyPress(1:lastTrial,1);
numberOfMissed = length(find(keyPress(1:lastTrial,1)==0));
allCorrect = find(response == presented);
numberOfCorrect = sum(response == presented)
performance = numberOfCorrect/lastTrial


%% performance of 1s trials
% load('M:\Data\phasecode\behavioral\pilot1.mat')

numberOfBlocks = 10;
trialsPerBlock = 80;
lastTrial = numberOfBlocks * trialsPerBlock;
idx1sec = find(conditions(1:lastTrial,6)==1);
presented = conditions(idx1sec,2);
response = keyPress(idx1sec,1);
numberOfMissed = length(find(keyPress(1:lastTrial,1)==0));
allCorrect = find(response == presented);
numberOfCorrect = sum(response == presented)

performance = numberOfCorrect/length(idx1sec)

%% performance over time
k=1;
Tperformance = zeros(20,1);
for b = 40:40:800
    Tidx1sec = find(conditions(b-39:b,6)==1);
    Tpresented = conditions(Tidx1sec,2);
    Tresponse = keyPress(Tidx1sec,1);
    numberOfCorrect = sum(Tresponse == Tpresented);
    Tperformance(k) = numberOfCorrect/40;
    k=k+1;
end
plot(Tperformance)

%% take the same trials as for the MEG data.
% 
% cd('/home/electromag/matves/Data/phasecode/behavioral');
% load('pilot1.mat');
% cd '/home/electromag/matves/MATLAB'
% numberOfBlocks = 4;
% 
% % cd('/home/electromag/matves/Data/phasecode/behavioral');
% % load('pilot_mats.mat');
% % cd '/home/electromag/matves/MATLAB'
% % numberOfBlocks = 2;
% 
% 
% behavior = [conditions keyPress]; % put conditions and keyPress together
% 
% % select presented conditions (if experiment is aborted before end, last
% % conditions trial will contain no response)
% trialsPerBlock = 80;
% lastTrial = numberOfBlocks * trialsPerBlock;
% behavior = behavior(1:lastTrial,:);
% 
% 
% % remove first trial in behavior. This trial wasn't saved in meg. ONLY FOR
% % PILOT1
% behavior = behavior(2:end, :);
% 
% % select only trials of 1sec
% idx1sec = find(behavior(:,6)==1);
% behavior = behavior(idx1sec,:);
% 
% % add the trial index
% behavior=[behavior zeros(length(behavior),1)];
% for j=1:length(behavior)
% behavior(j,9) = j;
% end
% 
% behavCorrect = behavior(behavior(:,2)==behavior(:,7), :);
% 
% 
% behavCorrect(4,:) = [];% too slow in meg.
% %% select correct meg trials
% meg = cfg.trl;
% meg = meg(:,4:end);
% megCorrect = meg(meg(:,7)==1, :);
% 
% 
% %% plot
% 
% %plot cues
% figure;
% subplot(3,1,1)
% bar(behavior(:,1))
% subplot(3,1,2)
% bar(meg(:, 1))
% subplot(3,1,3)
% bar(meg(:,1)-behavior(:,1))
% 
% %plot rotation
% figure;
% subplot(3,1,1)
% bar(behavior(:,2))
% subplot(3,1,2)
% bar(meg(:,3))
% subplot(3,1,3)
% bar(meg(:,3)-behavior(:,2))
% 
% 
% %plot response
% figure;
% subplot(3,1,1)
% bar(behavior(:,7))
% subplot(3,1,2)
% bar(meg(:, 5))
% subplot(3,1,3)
% bar(meg(:,5)-behavior(:,7))
% 
