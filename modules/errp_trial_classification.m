clearvars; clc;

subject = 'a5';
includepat  = {subject};
%includepat  = {'errp'};
excludepat  = {};
spatialfilter = 'car';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];
navpath     = 'analysis/navigation/';

datafiles = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
ndatafiles = length(datafiles);
util_bdisp(['[io] - Found ' num2str(ndatafiles) ' data files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

navfiles = util_getfile3(navpath, '.mat', 'include', includepat, 'exclude', excludepat);
nnavfiles = length(navfiles);
util_bdisp(['[io] - Found ' num2str(nnavfiles) ' navigation files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);


%% Data parameters
CommandTyp = 123;
ErrorTyp   = 987;

StartTyp = 100;
StopTyp  = 200;

%% Import data
util_bdisp(['[io] - Importing ' num2str(ndatafiles) ' files from ' datapath ':']);
[P, events, labels, classifiers, settings] = errp_concatenate_bandpass(datafiles);
nsamples  = size(P, 1);
nchannels = size(P, 2);
SampleRate = settings.data.samplerate;

%% Import navigation data
util_bdisp(['[io] - Importing ' num2str(nnavfiles) ' files from ' navpath ':']);
[cmd, pos, vel, navevents, navlabels] = errp_concatenate_navigation(navfiles);
nnavsamples = size(cmd, 1);

%% Filtering vel data
vz = errp_lowpass_velocity(vel(:, 3), SampleRate, 0.5);
cmdz = cmd(:, 3);
[joyk, evtJoy] = proc_get_event2([123 987], nnavsamples, navevents.POS, navevents.TYP, 1);
[vzthidx, joyidx] = find_cmd_latency(vz, joyk, 0.05);
%latency = vzthidx - joyidx;
latency = 0;

%% Downsampling
do_downsample = true;
SelectedSamplingRate = 64;
DownFactor = settings.data.samplerate/SelectedSamplingRate;
if do_downsample == true
    F = nan(size(P, 1)/DownFactor, size(P, 2));
    for chId = 1:size(P, 2)
        F(:, chId) = decimate(P(:, chId), DownFactor, 'FIR');
    end
    sRk = labels.samples.Rk(1:DownFactor:size(P, 1));
    events.POS = floor(events.POS/DownFactor);
    events.DUR = floor(events.DUR/DownFactor);
    SampleRate = SelectedSamplingRate;
else
    F = P;
    SampleRate = settings.data.samplerate;
end

%% Extract labels and events
%events = fix_start_duration(o_events, StartTyp, StopTyp);
%[RacK, evtRac] = proc_get_event2(StartTyp, nsamples, events.POS, events.TYP, events.DUR);
[CmdK, evtCmd] = proc_get_event2([CommandTyp ErrorTyp], nsamples, events.POS, events.TYP, 1);


%% Trial extraction
%latency = 0;
TrialPeriod = [0 1];
t_period = TrialPeriod(1):1/SampleRate:TrialPeriod(2) - 1/SampleRate;
ntrials = length(evtCmd.POS);
StartPos = latency + evtCmd.POS - floor(abs(TrialPeriod(1)*SampleRate));
StopPos  = latency + evtCmd.POS + floor(abs(TrialPeriod(2)*SampleRate)) - 1;
nsamples = length(StartPos(1):StopPos(1));
Trials = zeros(nsamples, nchannels, ntrials);
Sk = nan(ntrials, 1);
Rk = nan(ntrials, 1);
for trId = 1:ntrials
    cstart = StartPos(trId);
    cstop  = StopPos(trId);
    Trials(:, :, trId) = F(cstart:cstop, :);
    Sk(trId) = unique(labels.samples.Sk(cstart:cstop));
    Rk(trId) = unique(sRk(cstart:cstop));
end
TrialTYP = ones(ntrials, 1);
TrialTYP(evtCmd.TYP == ErrorTyp) = 0;

%% Temporal Fisher Score
K = nan(ntrials, nsamples*nchannels);

% Reshaping dataset for fisher score
for trId = 1:ntrials
    K(trId, :) = reshape(Trials(:, :, trId), nsamples*nchannels, 1);
end

% Computing fisher score
fs_t = proc_fisher2(K, TrialTYP);

%% Features extraction and classification on the whole dataset
WinLength = 10;
WinStart = 1:WinLength:nsamples-WinLength;
WinStop  = WinStart + WinLength;
nwindows = length(WinStart);

nfeatures = 10;
rW = extract_temporal_features(Trials, [WinStart; WinStop]);
[fs, fs_sel_id, fs_sel_val] = select_features(rW, TrialTYP, nfeatures);

[seltime, selchan] = proc_bin2bc([nwindows nchannels], fs_sel_id);

Model = fitcdiscr(rW(:, fs_sel_id), TrialTYP, 'DiscrimType','linear');
    
[Gk, pp] = predict(Model, rW(:, fs_sel_id));

TotAccuracy(1) = sum(Gk == TrialTYP)./ntrials;
TotAccuracy(2) = sum(Gk(TrialTYP == 1) == TrialTYP(TrialTYP == 1))./sum(TrialTYP == 1);
TotAccuracy(3) = sum(Gk(TrialTYP == 0) == TrialTYP(TrialTYP == 0))./sum(TrialTYP == 0);

%% One-leave out (run based)
runs = unique(Rk);
nruns = max(Rk);
TrainAccuracy = nan(3, nruns);
TestAccuracy = nan(3, nruns);
x_tr = cell(nruns, 1);
y_tr = cell(nruns, 1);
x_te = cell(nruns, 1);
y_te = cell(nruns, 1);
auc_tr = nan(nruns, 1);
auc_te = nan(nruns, 1);

for rId = 1:nruns
    cindex = Rk == runs(rId);

    % Train and test index
    trainIdx = ~cindex;
    testIdx  = cindex;

    rTrain = extract_temporal_features(Trials(:, :, trainIdx), [WinStart; WinStop]);
    rTest  = extract_temporal_features(Trials(:, :, testIdx), [WinStart; WinStop]);
    
    kTrain = TrialTYP(trainIdx);
    kTest  = TrialTYP(testIdx);

    % Select features and train model
    [fs, fs_sel_id, fs_sel_val] = select_features(rTrain, kTrain, nfeatures);
    Model = fitcdiscr(rTrain(:, fs_sel_id), kTrain, 'DiscrimType','quadratic');
    
    % Train accuracy
    [GkTr, ppTr] = predict(Model, rTrain(:, fs_sel_id));
    ntrials_tr = sum(trainIdx);
    TrainAccuracy(1, rId) = sum(GkTr == kTrain)./ntrials_tr;
    TrainAccuracy(2, rId) = sum(GkTr(kTrain == 1) == kTrain(kTrain == 1))./sum(kTrain == 1);
    TrainAccuracy(3, rId) = sum(GkTr(kTrain == 0) == kTrain(kTrain == 0))./sum(kTrain == 0);

    [x, y, t, auc] = perfcurve(kTrain, ppTr(:, 1), 0);

    x_tr{rId} = x;
    y_tr{rId} = y;
    auc_tr(rId) = auc;

    % Test Accuracy
    [GkTe, ppTe] = predict(Model, rTest(:, fs_sel_id));
    ntrials_te = sum(testIdx);
    TestAccuracy(1, rId) = sum(GkTe == kTest)./ntrials_te;
    TestAccuracy(2, rId) = sum(GkTe(kTest == 1) == kTest(kTest == 1))./sum(kTest == 1);
    TestAccuracy(3, rId) = sum(GkTe(kTest == 0) == kTest(kTest == 0))./sum(kTest == 0);

    [x, y, t, auc] = perfcurve(kTest, ppTe(:, 1), 0);
    x_te{rId} = x;
    y_te{rId} = y;
    auc_te(rId) = auc;
end


%% Figures
fig1 = figure;
fig_set_position(fig1, 'All');

layout     = 'eeg.antneuro.32.noeog.mi';
[~, ChannelList] = proc_get_montage(layout);

nrows = 3;
ncols = 2;

subplot(nrows, ncols, 1);
imagesc(reshape(fs_t, nsamples, nchannels)');
xtick = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(t_period(xtick)', '%4.2f'));
set(gca, 'YTick', 1:nchannels);
set(gca, 'YTickLabel', ChannelList);
xlabel('[s]');
ylabel('channel');
title('fisher score over time');
colorbar;

subplot(nrows, ncols, 2);

imagesc(reshape(fs, nwindows, nchannels)');
hold on;
plot(seltime, selchan, 'ko', 'MarkerFaceColor','k', 'MarkerSize',4);
hold off;
set(gca, 'YTick', 1:nchannels);
set(gca, 'YTickLabel', ChannelList);
xlabel('[window]');
ylabel('channel');
title('fisher score over window');
colorbar;

subplot(nrows, ncols, 3);
bar(100*TotAccuracy);
grid on;
ylim([0 100]);
set(gca, 'XTickLabel', {'Total', 'Correct', 'Error'});
ylabel('[%]');
title('Accuracy on the whole dataset');

subplot(nrows, ncols, 4);
hold on;
bar(100*[mean(TrainAccuracy, 2) mean(TestAccuracy, 2)]);
errorbar([(1:3) - 0.14; (1:3) + 0.14]', 100*[mean(TrainAccuracy, 2) mean(TestAccuracy, 2)], [], 100*[std(TrainAccuracy, [], 2) std(TestAccuracy, [], 2)], '.k');
hold off;
grid on;
ylim([0 100]);
set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', {'Total', 'Correct', 'Error'});
legend({'trainset', 'testset'});
ylabel('[%]');
title('Accuracy on the leave-one-out (per run)');

subplot(nrows, ncols, 5);
hold on;
for rId = 1:nruns
    plot(x_tr{rId}, y_tr{rId});
end
hold off;
grid on;
axis square;
title(['ROC on trainset - AUC: ' num2str(mean(auc_tr), '%4.2f') '+/-' num2str(std(auc_tr),'%4.2f')]);


subplot(nrows, ncols, 6);
hold on;
for rId = 1:nruns
    plot(x_te{rId}, y_te{rId});
end
hold off;
grid on;
axis square;
title(['ROC on testset - AUC: ' num2str(mean(auc_te), '%4.2f') '+/-' num2str(std(auc_te),'%4.2f')]);

sgtitle(['subject ' subject])

%% Functions
function rW = extract_temporal_features(P, windows)

    nwindows  = size(windows, 2);
    nchannels = size(P, 2);
    ntrials   = size(P, 3);

    W = nan(nwindows, nchannels, ntrials);

    for wId = 1:nwindows
        cstart = windows(1, wId);
        cstop  = windows(2, wId);
        W(wId, :, :) = mean(P(cstart:cstop, :, :));
    end
   
    rW = nan(ntrials, nwindows*nchannels);
    for trId = 1:ntrials
        rW(trId, :) = reshape(W(:, :, trId), nwindows*nchannels, 1);
    end
end

function [fs, fs_sel_id, fs_sel_val] = select_features(W, Wk, nfeatures)
    
    fs = proc_fisher2(W, Wk);

    [~, fs_w_sorted_idx] = sort(fs, 'descend');


    fs_sel_id = fs_w_sorted_idx(1:nfeatures);
    fs_sel_val = fs(fs_sel_id);
    
end