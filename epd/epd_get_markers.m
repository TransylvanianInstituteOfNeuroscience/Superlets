function markers = epd_get_markers(dataset, event_group_string)

event_group = parse_event_group_string(event_group_string);
markers     = [];

if (isfield(dataset, 'trial_count'))
    
    % go trial by trial searching for the event group match
    for i = 1 : dataset.trial_count
        [good, m] = get_event_group(dataset.trials(i).markers, event_group);
        if (good)
            markers = [markers, m];
        end
    end
    
else
    if (isfield(dataset, 'event_count'))
        
        % load the markers from file
        epd_markers = epd_load_events(dataset);
        
        % go through the events and validate them
        for i = 1 : numel(epd_markers)
            if (has_code(event_group, epd_markers(i).code))
                markers = [markers, epd_markers(i)];
            end
        end
    else
        error('input is neither an EPD or a trial structure');
    end
end
return

% parse a string describing an event group
function event_group = parse_event_group_string(event_group_string)

group       = epd_parse_event_pattern(event_group_string);
event_group = group{1};

return;

% check if an event group has a specific event code
function ok = has_code(event_group, code)

    ok = false;
    for i = 1 : numel(event_group)
        if (event_group(i) == code)
            ok = true;
        end
    end
    
return;

% match a trial marker array with an event group
function [ok, marker] = get_event_group(trial_markers, event_group)

ok      = false;
marker  = [];

for i = 1 : numel(trial_markers)
    for j = 1 : numel(event_group)
        if (trial_markers(i).code == event_group(j))
            ok      = true;
            marker  = trial_markers(i);
            return;
        end
    end
end
return;


