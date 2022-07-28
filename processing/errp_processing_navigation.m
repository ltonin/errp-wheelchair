clearvars; clc;

subject = 'c5';

includepat  = {subject};
excludepat  = {};
depthlevel  = 3;

rootpath    = '/mnt/data/Research/';
folder      = 'errp_wheelchair';
gdfpath     = [rootpath '/' folder '/'];
savedir     = 'analysis/navigation/';
recompute   = true;

SampleRate = 512;

%% Get datafiles
files = util_getfile3(gdfpath, '.bag', 'include', includepat, 'exclude', excludepat, 'level', depthlevel);

nfiles = length(files);
if(nfiles > 0)
    util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
else
    error(['[io] - No files found with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
end

%% Create/Check for savepath
util_mkdir(pwd, savedir);

%% Processing files

for fId = 1:nfiles

    cfullname = files{fId};
    [cfilepath, cfilename, cfileext] = fileparts(cfullname);
    
    util_bdisp(['[io] + Loading file ' num2str(fId) '/' num2str(nfiles)]);
    disp(['     |-File: ' cfullname]);
    
    %% Check if the file has been already processed
    [~, pfilename] = fileparts(cfullname);
    if (recompute == false) && exist([savedir pfilename '.mat'], 'file') == 2
        disp('     |-Processed bandpass already exists. Skipping the recomputing');
        continue;
    end

    %% Loading file
    cbag = rosbag(cfullname);

    sel_cmdvel = select(cbag, 'Topic', '/cmd_vel');
    msg_cmdvel = readMessages(sel_cmdvel, 'DataFormat', 'struct');
    
    sel_evtbus = select(cbag, 'Topic', '/events/bus');
    msg_evtbus = readMessages(sel_evtbus, 'DataFormat', 'struct');

    sel_joy = select(cbag, 'Topic', '/joy');
    msg_joy = readMessages(sel_joy, 'DataFormat', 'struct');
    
    sel_odom = select(cbag, 'Topic', '/odom');
    msg_odom = readMessages(sel_odom, 'DataFormat', 'struct');

    %% Extracting messages
    % X and Z command velocities
    X_cmdvel = cellfun(@(m) double(m.Linear.X), msg_cmdvel);
    Z_cmdvel = cellfun(@(m) double(m.Angular.Z), msg_cmdvel);
    t_cmdvel = sel_cmdvel.MessageList.Time - cbag.StartTime;

    % X and Z positions and odom velocities
    pX_odom   = cellfun(@(m) double(m.Pose.Pose.Position.X), msg_odom);
    pY_odom   = cellfun(@(m) double(m.Pose.Pose.Position.Y), msg_odom);
    pQ_odom   = cellfun(@(m) double(quat2eul(rosReadQuaternion(m.Pose.Pose.Orientation))), msg_odom, 'UniformOutput',false);
    pQ_odom   = cell2mat(pQ_odom);
    pZ_odom   = pQ_odom(:, 1);
    vX_odom   = cellfun(@(m) double(m.Twist.Twist.Linear.X), msg_odom);
    vZ_odom   = cellfun(@(m) double(m.Twist.Twist.Angular.Z), msg_odom);
    t_odom   = sel_odom.MessageList.Time - cbag.StartTime;

    % Joy events
    buttons_joy_cell = cellfun(@(m) double(m.Buttons), msg_joy, 'UniformOutput',false);
    buttons_joy = [];
    t_joy = [];
    for bId = 1:length(buttons_joy_cell)
        currIdx = find(buttons_joy_cell{bId});

        if(isempty(currIdx) == false)
            buttons_joy = cat(1, buttons_joy, currIdx);
            t_joy = cat(1, t_joy, sel_joy.MessageList.Time(bId)*ones(length(currIdx), 1));
        end
    end
    t_joy = t_joy - cbag.StartTime;



    % Event bus
    id_evt  = cellfun(@(m) double(m.Event), msg_evtbus);
    t_evt   = sel_evtbus.MessageList.Time - cbag.StartTime;

    %% Converting analog signals in timeseries
    ts_cX = timeseries(X_cmdvel, t_cmdvel);
    ts_cZ = timeseries(Z_cmdvel, t_cmdvel);

    ts_pX = timeseries(pX_odom, t_odom);
    ts_pY = timeseries(pY_odom, t_odom);
    ts_pZ = timeseries(pZ_odom, t_odom);
    ts_vX = timeseries(vX_odom, t_odom);
    ts_vZ = timeseries(vZ_odom, t_odom);

    t = 0:1/SampleRate:ceil(cbag.EndTime - cbag.StartTime);

    %% Resample analog signals with respect to uniform time support
    tts_cX = resample(ts_cX, t, 'zoh');
    tts_cZ = resample(ts_cZ, t, 'zoh');
    tts_pX = resample(ts_pX, t, 'zoh');
    tts_pY = resample(ts_pY, t, 'zoh');
    tts_pZ = resample(ts_pZ, t, 'zoh');
    tts_vX = resample(ts_vX, t, 'zoh');
    tts_vZ = resample(ts_vZ, t, 'zoh');
    
    %% Create events structure and synchronized with time support
    events_bus = create_event_structure(id_evt, t_evt, t, SampleRate);
    events_joy = create_event_structure(buttons_joy, t_joy, t, SampleRate);

    events = merge_event_structure(events_bus, events_joy);

    %% Renaming analog variables
    cmd_x = tts_cX.Data;
    cmd_z = tts_cZ.Data;
    pos_x = tts_pX.Data;
    pos_y = tts_pY.Data;
    pos_z = tts_pZ.Data;
    vel_x = tts_vX.Data;
    vel_z = tts_vZ.Data;
    

    %% Saving the navigation data
    sfilename = [savedir '/' pfilename '.mat'];
    util_bdisp(['[out] - Saving navigation data in: ' sfilename]);
    save(sfilename, 'cmd_x', 'cmd_z', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_z', 't', 'events', 'SampleRate'); 


end

function events = create_event_structure(eTYP, eTIM, t_support, samplerate)

    nevt = length(eTYP);
    if(nevt ~= length(eTIM))
        error('chk:len: TYP and TIM must have the same dimension')
    end

    if (min(eTIM) < min(t_support) || max(eTIM) > max(t_support))
        error('chk:tim: TIM provided is not within t_support')
    end

    TYP = nan(nevt, 1);
    DUR = nan(nevt, 1);
    POS = nan(nevt, 1);
    ERR = nan(nevt, 1);

    for eId = 1:nevt
        cevt_time = eTIM(eId);
        cevt_typ  = eTYP(eId);

        [cevt_err, cevt_err_idx] = min(abs(cevt_time - t_support));

        TYP(eId) = cevt_typ;
        POS(eId) = cevt_err_idx;
        DUR(eId) = 0;
        ERR(eId) = (cevt_err*samplerate);

    end

    events.TYP = TYP;
    events.POS = POS;
    events.DUR = DUR;
    events.ERR = ERR;

end

function events = merge_event_structure(evt1, evt2) 

    TYP1 = evt1.TYP;
    POS1 = evt1.POS;
    DUR1 = evt1.DUR;
    ERR1 = evt1.ERR;

    TYP2 = evt2.TYP;
    POS2 = evt2.POS;
    DUR2 = evt2.DUR;
    ERR2 = evt2.ERR;


    TYP = [TYP1; TYP2];
    DUR = [DUR1; DUR2];
    ERR = [ERR1; ERR2];

    [POS, idx] = sort([POS1; POS2]);
    TYP = TYP(idx);
    DUR = DUR(idx);
    ERR = ERR(idx);

    events.TYP = TYP;
    events.POS = POS;
    events.DUR = DUR;
    events.ERR = ERR;
    


end