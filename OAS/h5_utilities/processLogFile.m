function [log_events] = processLogFile(foldername,time, isremote)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

mmPerStep = 0.001253814; % calibration for gcamp + behavior tracker

%% process log file
[fld, ~, ~]=fileparts(foldername);
logd = dir([fld '\*.txt']);
acq = 0; % acquires log data within recording range when set to 1;
stimTimes = [];
xLoc = NaN(length(time),1);
yLoc = NaN(length(time),1);
starttime = time(1);

for i = 1:length(logd)
    log_path = fullfile(logd(i).folder, logd(i).name);

    if isremote == 1
        local_path = ['C:\tmp\' logd(i).name];
        copyfile(log_path, local_path)
        fid = fopen(local_path,"r");
    elseif isremote == 0
        fid = fopen(log_path,"r");
    end

    while~feof(fid)
        line = fgetl(fid);
        l = regexp(line, ' ', 'split');
        lTime = str2double(l{1}); % time at current line

        if contains(line, 'command sent: DO _writer_start') % find a line matching the beginning of the recording
            ls = regexp(line, ' ', 'split');
            logStart = str2double(ls{1});
            if abs(logStart-starttime)<1
                acq= 1;
                disp(line)
            end
        end

        % stop acquisition when recording ends
        if lTime>max(time)
            disp(line)
            acq = 0;
            break
        end

        if acq == 1
            if contains(line, '<CLIENT WITH GUI> command sent: DO _teensy_commands_set_led o 1')
                stimTimes(end+1,1) = alignEvent(line,time);
                disp(line)
            end
        end

        if acq == 1
            if contains(line, 'tracker_behavior received position')
                locTime = alignEvent(line,time);
                r = regexp(line,' ', 'split');
                r = regexp(r{end}, '(-?\d+)', 'match');
                xl = str2double(r{1})*mmPerStep; % X coordinate in mm units
                yl = str2double(r{2})*mmPerStep; % Y coordinate in mm units

                xLoc(locTime,1) = xl;
                yLoc(locTime,1) = yl;




            end
        end

    end
    fclose(fid);

    if isremote == 1
        delete(local_path)
    end

end

%% calculate instantaneous velocity

reltime = time-starttime;
firstsec = find(reltime>1,1);
secondsec = find(reltime>2,1);

sec = secondsec-firstsec;
velocity =NaN(length(time),1);

for i = 2:length(xLoc)-(sec+1)
    dx = xLoc(i)-xLoc(i+sec); %change in xLoc per second
    dy = yLoc(i)-yLoc(i+sec); %change in yLoc per second
    velocity(i) = sqrt(dx.^2 + dy.^2);
end

log_events.time = time;
log_events.startime = starttime;
log_events.stimTimes =stimTimes;
log_events.velocity = velocity;
log_events.xLoc = xLoc;
log_events.yLoc = yLoc;

end


function [idx] = alignEvent(event, time)
et = regexp(event, ' ', 'split');
eTime = str2double(et{1});
idx = find(time>=eTime,1);
end