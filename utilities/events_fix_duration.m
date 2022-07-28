function events = events_fix_duration(events, PadEvents, RaceEvents)
% EVENTS_FIX_DURATION fixes the duration of the Pad and Race events

    POS = events.POS;
    TYP = events.TYP;
    DUR = events.DUR;

    % Fix Pad duration
    RaceStart = find(TYP == RaceEvents(1));
    RaceStop  = find(TYP == RaceEvents(2));

    cRk = zeros(length(POS), 1);
    for rId = 1:length(RaceStart)
        cstart = RaceStart(rId);
        cstop  = RaceStop(find(RaceStop > cstart, 1));
        cRk(cstart:cstop) = rId;
    end

    npads = length(PadEvents);
    index = false(length(TYP), 1);
    for pId = 1:npads
        index = index | TYP == PadEvents(pId);
    end

    for rId = 1:length(RaceStart)
        evtIdx = find(index & cRk == rId);
        DUR(evtIdx) = [POS(evtIdx(2:end)); POS(find(cRk == rId, 1, 'last'))] - POS(evtIdx) - 1;
    end

    % Fix Race event

    real_race_start = diff(cRk > 0) == 1;
    real_race_stop  = diff(cRk > 0) == -1;
    DUR(TYP == RaceEvents(1)) = POS(real_race_stop) - POS(real_race_start) - 1;
        
    incIdx = setdiff(1:length(POS), find(TYP == RaceEvents(2)));
    POS = POS(incIdx);
    DUR = DUR(incIdx);
    TYP = TYP(incIdx);
    
    events.POS = POS;
    events.DUR = DUR;
    events.TYP = TYP;
end