clearvars; clc;

subject = 'c5';
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
F = P;

%% Import navigation data
util_bdisp(['[io] - Importing ' num2str(nnavfiles) ' files from ' navpath ':']);
[cmd, pos, vel, navevents, navlabels] = errp_concatenate_navigation(navfiles);
nnavsamples = size(cmd, 1);

%% Filtering vel data
vz = errp_lowpass_velocity(vel(:, 3), SampleRate, 0.5);
cmdz = cmd(:, 3);


[joyk, evtJoy] = proc_get_event2([123 987], nnavsamples, navevents.POS, navevents.TYP, 1);

[vzthidx, joyidx] = find_cmd_latency(vz, joyk, 0.05);
latency = vzthidx - joyidx;
latency = 0;
%% Downsampling
% do_downsample = true;
% SelectedSamplingRate = 64;
% DownFactor = settings.data.samplerate/SelectedSamplingRate;
% if do_downsample == true
%     F = nan(size(P, 1)/DownFactor, size(P, 2));
%     for chId = 1:size(P, 2)
%         F(:, chId) = decimate(P(:, chId), DownFactor, 'FIR');
%     end
%     events.POS = floor(events.POS/DownFactor);
%     events.DUR = floor(events.DUR/DownFactor);
%     SampleRate = SelectedSamplingRate;
% else
%     F = P;
%     SampleRate = settings.data.samplerate;
% end

%% Extract labels and events
%events = fix_start_duration(o_events, StartTyp, StopTyp);
%[RacK, evtRac] = proc_get_event2(StartTyp, nsamples, events.POS, events.TYP, events.DUR);
[CmdK, evtCmd] = proc_get_event2([CommandTyp ErrorTyp], nsamples, events.POS, events.TYP, 1);

% if strcmpi(subject, 'a7')
%     Offset = -1.5;
%     evtCmd.POS = evtCmd.POS + floor(Offset*SampleRate);
% end


%% Trial extraction
%latency = 0;
TrialPeriod = [-0.5 1.0];
ntrials = length(evtCmd.POS);
StartPos = latency + evtCmd.POS - floor(abs(TrialPeriod(1)*SampleRate));
StopPos  = latency + evtCmd.POS + floor(abs(TrialPeriod(2)*SampleRate)) - 1;
Trials = zeros(length(StartPos(1):StopPos(1)), nchannels, ntrials);
Sk = nan(ntrials, 1);
for trId = 1:ntrials
    cstart = StartPos(trId);
    cstop  = StopPos(trId);
    Trials(:, :, trId) = F(cstart:cstop, :);
    Sk(trId) = unique(labels.samples.Sk(cstart:cstop));
end
TrialTYP = ones(ntrials, 1);
TrialTYP(evtCmd.TYP == ErrorTyp) = 0;

%% Figures
fig1 = figure;
fig_set_position(fig1, 'All');

layout     = 'eeg.antneuro.32.noeog.mi';
[~, ChannelList] = proc_get_montage(layout);
%ChannelLbs = {'Fz', 'FCz', 'Cz', 'Fp1', 'Fp2'};
ChannelLbs = {'Fz', 'FCz', 'Cz'};
ChannelIds = proc_get_channel(ChannelLbs, ChannelList);
NumChannels = length(ChannelIds);

AvgCorrect = mean(Trials(:, :, TrialTYP == 1), 3);
AvgError   = mean(Trials(:, :, TrialTYP == 0), 3);

t = TrialPeriod(1):1/SampleRate:TrialPeriod(2) - 1/SampleRate;

hchans = zeros(NumChannels, 1);
for chId = 1:NumChannels
    subplot(3, NumChannels, chId);

    hold on;
    plot(t, AvgCorrect(:, ChannelIds(chId)), 'b', 'LineWidth',2);
    plot(t, AvgError(:, ChannelIds(chId)), 'r', 'LineWidth',2);
    plot(t, AvgError(:, ChannelIds(chId)) - AvgCorrect(:, ChannelIds(chId)), 'k', 'LineWidth',2);
    hold on;
    %ylim([-6 6]);
    grid on;
    plot_vline(0, 'k');
    ylabel('amplitude [uV]');
    xlabel('time [s]')
    title(ChannelLbs{chId});
    legend('Correct', 'Error', 'Error - Correct', 'Location','best');
    hchans(chId) = gca;
end
plot_set_limits(hchans, 'y', 'minmax');

himge = nan(NumChannels, 1);
for chId = 1:NumChannels
    subplot(3, NumChannels, NumChannels + chId);
    imagesc(t, 1:sum(TrialTYP == 0), squeeze(Trials(:, ChannelIds(chId), TrialTYP == 0))');
    plot_vline(0, 'w')
    %colormap('bone');
    %caxis([-10 10]);
    title('Error')
    himge(chId) = gca;
end

himgc = nan(NumChannels, 1);
for chId = 1:NumChannels
    subplot(3, NumChannels, 2*NumChannels + chId);
    imagesc(t, 1:sum(TrialTYP == 1), squeeze(Trials(:, ChannelIds(chId), TrialTYP == 1))');
    plot_vline(0, 'w')
    %colormap('bone');
    %caxis([-10 10]);
    title('Correct')
    himgc(chId) = gca;
end

himg = [himge; himgc];
cmin = min(min(cell2mat(get(himg, 'CLim'))));
cmax = max(max(cell2mat(get(himg, 'CLim'))));

set(himg, 'CLim', unique(cell2mat(get(hchans, 'YLim')))*5);
set(himg, 'CLim', [-10 10])

%% Topoplots
fig2 = figure;
fig_set_position(fig2, 'Top');

layout     = 'eeg.antneuro.32.noeog.mi';
[~, ChannelList] = proc_get_montage(layout);

NumChannels = length(ChannelIds);


TopoPeriods(:, 1) = [-0.05 0.05];
TopoPeriods(:, 2) = [0.2 0.3];
TopoPeriods(:, 3) = [0.5 0.7];
TopoPeriods(:, 4) = [0.7 1.0];

nperiods = size(TopoPeriods, 2);
load('chanlocs64.mat');
% newchanlocs = convert_chanlocs(chanlocs, ChannelList);
newchanlocs = chanlocs;

for tId = 1:nperiods
    [~, cstart] = min(abs(t - TopoPeriods(1, tId)));
    [~, cstop] = min(abs(t - TopoPeriods(2, tId)));
    
    subplot(2, nperiods, tId);
    [cdata_err, cindex_err] = errp_fill_layout(mean(AvgError(cstart:cstop, :)), ChannelList, newchanlocs);
    topoplot(cdata_err, newchanlocs, 'maplimits', unique(cell2mat(get(hchans, 'YLim'))), 'headrad', 'rim', 'conv', 'on', 'shading', 'interp', 'electrodes', 'off', 'emarker2', {cindex_err, '.', 'k'});
    title([num2str(TopoPeriods(1, tId)) '-' num2str(TopoPeriods(2, tId)) ' s']);
    axis image
    subplot(2, nperiods, tId + nperiods);
    [cdata_corr, cindex_corr] = errp_fill_layout(mean(AvgCorrect(cstart:cstop, :)), ChannelList, newchanlocs);
    topoplot(cdata_corr, newchanlocs, 'maplimits', unique(cell2mat(get(hchans, 'YLim'))), 'headrad', 'rim', 'conv', 'on', 'shading', 'interp', 'electrodes', 'off', 'emarker2', {cindex_corr, '.', 'k'});
    title([num2str(TopoPeriods(1, tId)) '-' num2str(TopoPeriods(2, tId)) ' s']);
    axis image
end



%% Functions
function events = fix_start_duration(events, StartEvt, StopEvt)
    
    events.DUR(events.TYP == StartEvt) = events.POS(events.TYP == StopEvt) - events.POS(events.TYP == StartEvt) + 1;

end

function nchanlocs = convert_chanlocs(chanlocs, chanlist)

    %nchanlocs = nan(length(chanlist), 1);
    for lId = 1:length(chanlist)
        cname = char(chanlist(lId));
        
        for chId = 1:length(chanlocs)
            if(strcmpi(chanlocs(chId).labels, cname) == true)
                nchanlocs(lId) = chanlocs(chId);
                break;
            end

            %disp(['cannot find channel ' cname ' in chanlocs'])
        end
    end

end

function [fdata, findex] = fill_layout(data, chanlist, locations)

    fdata = zeros(length(locations), 1);
    nchannels = size(data,2);

    chanindex = nan(nchannels, 1);
    for chId = 1:nchannels
        cname = char(chanlist(chId));

        for lId = 1:length(locations)
            if(strcmpi(locations(lId).labels, cname) == true)
                chanindex(chId) = lId;
                break;
            end
        end
    end
    
    fdata(chanindex) = data;
    findex = chanindex;

end
