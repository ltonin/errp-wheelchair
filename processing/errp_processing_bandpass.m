clearvars; clc;

subject = 'c5';

includepat  = {subject};
excludepat  = {};
depthlevel  = 2;

rootpath    = '/mnt/data/Research/';
folder      = 'errp_wheelchair';
gdfpath     = [rootpath '/' folder '/'];

artifactrej       = 'none'; % {'FORCe', 'none'}
ForceWinLength    = 1.0;
spatialfilter     = 'car';
savedir           = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];
recompute         = true;

%% Processing parameters
nchannels  = 32;
bands      = [1 10];
filtorder  = 4;
layout     = 'eeg.antneuro.32.noeog.mi';

%% Get datafiles
files = util_getfile3(gdfpath, '.gdf', 'include', includepat, 'exclude', excludepat, 'level', depthlevel);

NumFiles = length(files);
if(NumFiles > 0)
    util_bdisp(['[io] - Found ' num2str(NumFiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
else
    error(['[io] - No files found with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
end

%% Create/Check for savepath
util_mkdir(pwd, savedir);

%% Processing files
for fId = 1:NumFiles
    cfullname = files{fId};
    [cfilepath, cfilename, cfileext] = fileparts(cfullname);
    
    util_bdisp(['[io] + Loading file ' num2str(fId) '/' num2str(NumFiles)]);
    disp(['     |-File: ' cfullname]);
    
    %% Check if the file has been already processed
    [~, pfilename] = fileparts(cfullname);
    if (recompute == false) && exist([savedir pfilename '.mat'], 'file') == 2
        disp('     |-Processed bandpass already exists. Skipping the recomputing');
        continue;
    end
    
    %% Loading data
    disp('     |-Loading GDF data');
    try
        [s, h] = sload(cfullname);
        s = s(:, 1:nchannels);
    catch ME
        warning('[warning] - Cannot load filename. Skipping it.');
        warning(['[warning] - Error: ' ME.message]);
        continue;
    end
    
    
    %% Get information from filename
    cinfo = whi_util_get_info(cfullname);
    
    %% Get montage 
    [montage, labels] = proc_get_montage(layout);
    
    %% Processing data
    util_bdisp('[proc] + Processing the data');

    %s(:, 25) = (s(:, 9) + s(:, 10))./2;

    % DC
    disp('       |-Removing DC');
    s_dc = s - repmat(mean(s, 1), size(s, 1), 1);
    s_dc = s;

    % Compute Spatial filter
    disp(['       |-Spatial filter: ' spatialfilter ' (' layout ')']);

    switch(spatialfilter)
        case 'none'
            s_filt = s_dc;
        case 'car'
            s_filt = proc_car(s_dc);
        case 'laplacian'
            lap = whi_proc_laplacian_mask(montage, nchannels);
            s_filt = s_dc*lap;
        otherwise
            error(['Unknown spatial filter selected ' spatialfilter]);
    end
    
    % Compute bandpass filters
    s_bp = filt_bp(s_filt, filtorder, bands, h.SampleRate);

    P = s_bp;
    
    % Extracting events
    disp('       |-Extract events');
    cevents     = h.EVENT;
    events.TYP = cevents.TYP;
    events.POS = cevents.POS;
    events.DUR = cevents.DUR;
    
    % Modality
    disp('       |-Extract additional info (modality, protocol, date)');
    modality = cinfo.modality;
    
    % Protocol
    switch(cinfo.modality)
        case 'offline'
            protocol = 'bci-calibration';
        case 'online'
            protocol = 'bci-training';
        case {'feedback', 'control', 'navigation', 'cybathlon'}
            protocol = 'bci-race';
        otherwise
            protocol = 'unknown';
    end
    
    %% Create settings structure
    settings.data.filename          = cfullname;
    settings.data.nsamples          = size(s, 1);
    settings.data.nchannels         = size(s, 2);
    settings.data.lchannels         = labels;
    settings.data.samplerate        = h.SampleRate;
    settings.artifact.name          = artifactrej;
    settings.artifact.force.wlength = ForceWinLength;
    settings.spatial.filter         = spatialfilter;
    settings.bandpass.order         = filtorder;
    settings.bandpass.bands         = bands;
    settings.modality.legend        = {'offline','online', 'feedback', 'control', 'navigation', 'cybathlon'};
    settings.modality.name          = modality;
    settings.protocol.legend        = {'bci-calibration', 'bci-training', 'bci-race', 'unknown'};
    settings.protocol.name          = protocol;
    settings.info                   = cinfo;
    
    
    sfilename = [savedir '/' pfilename '.mat'];
    util_bdisp(['[out] - Saving bandpass in: ' sfilename]);
    save(sfilename, 'P', 'events', 'settings'); 
end

