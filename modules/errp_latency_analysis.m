clearvars; clc;

subject = 'c5';
includepat  = {subject};
excludepat  = {};
spatialfilter = 'car';
artifactrej   = 'none'; % {'FORCe', 'none'}
gdfpath    = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];
navpath    = 'analysis/navigation/';

gdffiles  = util_getfile3(gdfpath, '.mat', 'include', includepat, 'exclude', excludepat);
ngdffiles = length(gdffiles);
util_bdisp(['[io] - Found ' num2str(ngdffiles) ' GDF files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

navfiles  = util_getfile3(navpath, '.mat', 'include', includepat, 'exclude', excludepat);
nnavfiles = length(navfiles);
util_bdisp(['[io] - Found ' num2str(nnavfiles) ' navigation files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

%% General Parameters
CommandTyp = 123;
ErrorTyp   = 987;
StartTyp   = 100;
StopTyp    = 200;

SampleRate = 512;
NavLowPassCutOff = 0.5;

%% Import GDF data
util_bdisp(['[io] - Importing ' num2str(ngdffiles) ' GDF files from ' gdfpath ':']);
[~, gdfevents, gdflabels] = errp_concatenate_bandpass(gdffiles);
nsamplesGDF = length(gdflabels.samples.Rk);
runsGDF = unique(gdflabels.samples.Rk);
nrunsGDF = length(runsGDF);

%% Import navigation data
util_bdisp(['[io] - Importing ' num2str(nnavfiles) ' navigationfiles from ' navpath ':']);
[cmd, pos, vel, navevents, navlabels] = errp_concatenate_navigation(navfiles);
nsamplesNav = length(navlabels.Rk);
runsNav = unique(navlabels.Rk);
nrunsNav = length(runsNav);

%% Extracting GDF cmd events
[GDFCmdK, evtGDFCmd] = proc_get_event2([CommandTyp ErrorTyp], nsamplesGDF, gdfevents.POS, gdfevents.TYP, 1);

evtGDFCmd.RUN = nan(length(evtGDFCmd.POS), 1);
for rId = 1:nrunsGDF
    cindex = gdflabels.samples.Rk(evtGDFCmd.POS) == runsGDF(rId);
    evtGDFCmd.RUN(cindex) = runsGDF(rId);
end

%% Extracting Navigation events
[NavJoyK, evtNavJoy] = proc_get_event2([CommandTyp ErrorTyp 118], nsamplesNav, navevents.POS, navevents.TYP, 1);
evtNavJoy.RUN = nan(length(evtNavJoy.POS), 1);
for rId = 1:nrunsNav
    cindex = navlabels.Rk(evtNavJoy.POS) == runsNav(rId);
    evtNavJoy.RUN(cindex) = runsNav(rId);
end

%% Check equality between evtNav e evtGDF

if(isequal(evtNavJoy.TYP, evtGDFCmd.TYP) == false)
    error('GDF and navigation event TYP are different')
end

if(isequal(evtNavJoy.RUN, evtGDFCmd.RUN) == false)
    error('GDF and navigation event RUN are different')
end

runs = runsNav;
nruns = nrunsNav;

%% Extracting Navigation (valid) commands
cmdz_raw = [0; diff(cmd(:, 3) ~= 0) == 1];
njoy_raw = NavJoyK ~= 0;

njoy_raw_idx = find(NavJoyK);
njoy_padding = zeros(length(njoy_raw), 1);
for idx = 1:length(njoy_raw_idx)
    cstart = njoy_raw_idx(idx) - 5;
    cstop  = njoy_raw_idx(idx) + 5;
    njoy_padding(cstart:cstop) = 1;
end

cmdz = njoy_padding & cmdz_raw;
evtCmdz.POS = find(cmdz);
evtCmdz.TYP = evtNavJoy.TYP;
evtCmdz.DUR = 1;
evtCmdz.RUN = evtNavJoy.RUN;

%% Latency between CmdZ (valid) and JoyEvt per run

mlatency_cmdjoy = nan(nruns, 1);
slatency_cmdjoy = nan(nruns, 1);
for rId = 1:nruns
    cindex = evtCmdz.RUN == runs(rId);
    mlatency_cmdjoy(rId) = mean( (evtCmdz.POS(cindex) - evtNavJoy.POS(cindex))./SampleRate);
    slatency_cmdjoy(rId) = std( (evtCmdz.POS(cindex) - evtNavJoy.POS(cindex))./SampleRate);
end

%% Processing rotational velocity
velz_raw = vel(:, 3);
velz = errp_lowpass_velocity(velz_raw, SampleRate, NavLowPassCutOff);
VelThresholdPerc = 0.05;
latency_velz = nan(length(evtCmdz.POS), 1);
for eId = 1:length(evtCmdz.POS)
    cpos = evtCmdz.POS(eId);
    cstart = cpos;
    if eId == length(evtCmdz.POS)
        cstop = length(velz);
    else
        cstop  = evtCmdz.POS(eId+1);
    end
   
    cvel = velz(cpos);
    cthresholds = [cvel+cvel*VelThresholdPerc cvel-cvel*VelThresholdPerc];

    cvelIdx = find(velz(cstart:cstop) >= cthresholds(1) |  velz(cstart:cstop) <= cthresholds(2), 1, 'first');

    if(isempty(cvelIdx) == false)
        latency_velz(eId) = cpos + cvelIdx;
    end
end
mlatency_cmdvelz = nan(nruns, 1);
slatency_cmdvelz = nan(nruns, 1);
for rId = 1:nruns
    cindex = evtCmdz.RUN == runs(rId);
    cmdpos = evtCmdz.POS(cindex);
    velpos = latency_velz(cindex);
    mlatency_cmdvelz(rId) = mean(velpos-cmdpos);
    slatency_cmdvelz(rId) = std(velpos-cmdpos);
    
end

%% Extracting start event for GDF and Nav
NavStartK = proc_get_event2(StartTyp, nsamplesNav, navevents.POS, navevents.TYP, 1);
GdfStartK = proc_get_event2(StartTyp, nsamplesGDF, gdfevents.POS, gdfevents.TYP, 1);
latency_per_run = nan(nruns, 1);
for rId = 1:nruns
    cindex_nav = navlabels.Rk == runs(rId);
    cindex_gdf = gdflabels.samples.Rk == runs(rId);

    firstNav = find(NavStartK ~= 0 & cindex_nav, 1, 'first');
    firstGDF = find(GdfStartK ~= 0 & cindex_gdf, 1, 'first');
    latency_per_run(rId) = firstNav - firstGDF;
end

%% Figures
fig1 = figure;
fig_set_position(fig1, 'All');

nrows = 2;
ncols = 3;
run_index = [1; find(diff(navlabels.Rk) == 1)];

% Rotational velocity
subplot(nrows, ncols, [1 2 3]);
hold on;
plot(velz);
plot(evtCmdz.POS, velz(evtCmdz.POS), 'ko');
plot(latency_velz, velz(latency_velz), 'go');
set(gca, 'XTick', run_index);
set(gca, 'XTickLabel', num2str((1:nruns)'));
plot_vline(run_index, 'k');
grid on;
xlabel('time [run]');
ylabel('vel z [m/s]');

% Latency between CmdZ (valid) and JoyEvt per run
subplot(nrows, ncols, 4);
hold on;
bar(mlatency_cmdjoy);
errorbar(1:nruns, mlatency_cmdjoy, [], slatency_cmdjoy, 'k.');
lx = plot_hline(mean(mlatency_cmdjoy), 'r');
hold off;
set(gca, 'XTick', 1:nruns);
set(gca, 'XTickLabel', 1:nruns);
ylim([0 2/SampleRate]);
grid on;
ax = gca;
ax.YAxis.Exponent = 0;
ylabel('latency [s]');
xlabel('run');
title('Latency CmdVel vs. JoyEvt')
legend(lx, 'average')

% Latency between CmdZ (valid) and VelZ over threshold per run
subplot(nrows, ncols, 5);
hold on;
bar(mlatency_cmdvelz/SampleRate);
errorbar(1:nruns, mlatency_cmdvelz/SampleRate, [], slatency_cmdvelz/SampleRate, 'k.');
lx = plot_hline(mean(mlatency_cmdvelz/SampleRate), 'r');
hold off;
set(gca, 'XTick', 1:nruns);
set(gca, 'XTickLabel', 1:nruns);
grid on;
ax = gca;
ax.YAxis.Exponent = 0;
ylabel('latency [s]');
xlabel('run');
title(['Latency CmdVel vs. VelZ (Threshold: ' num2str(100*VelThresholdPerc) '%)']);
legend(lx, 'average')


