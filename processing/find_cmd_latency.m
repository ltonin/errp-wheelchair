function [vzThIndex, joyidx] = find_cmd_latency(vz, joyk, Threshold)


    joy = joyk ~= 0;

    joyidx = find(joy);
    
    vzAtjoy = vz(joyidx);

    vzThIndex = nan(length(joyidx), 1);
    for idx = 1:length(joyidx)
        try
            cstart = joyidx(idx);
            if(idx < length(joyidx))
                cstop  = joyidx(idx+1);
            else
                cstop = length(vz);
            end
    
            cvz = vz(cstart:cstop);
            
            cvzth_pos = cvz(1) + cvz(1)*Threshold;
            cvzth_neg = cvz(1) - cvz(1)*Threshold;
    
            cvzThIdx = find(cvz >= cvzth_pos | cvz <= cvzth_neg, 1, 'first');
            
            if isempty(cvzThIdx)
                keyboard
            else
                vzThIndex(idx) = cvzThIdx + cstart;
            end
        catch
            keyboard
        end

    end
    

end