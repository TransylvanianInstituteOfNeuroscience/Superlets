function trial_structure = epd_parse_trial_structure(epd, event_pattern_string, trial_info)

event_pattern               = epd_parse_event_pattern(event_pattern_string);
markers                     = epd_load_events(epd);
trial_markers               = [];
trial_structure.trials      = [];
trial_structure.trial_count = 0;

% parse the trial structure
for i = 1 : numel(markers)
    if (has_code(event_pattern{numel(trial_markers) + 1}, markers(i).code))
        trial_markers = [trial_markers, markers(i)];
    end
    
    if (numel(trial_markers) == numel(event_pattern))
        trial_structure.trial_count                                 = trial_structure.trial_count + 1;
        trial_structure.trials(trial_structure.trial_count).markers = trial_markers;
        trial_structure.trials(trial_structure.trial_count).info    = { };
        trial_markers                                               = [];
    end
end


% assign trial info to the trials
if (~isempty(trial_info))
    if (trial_info.trial_count ~= trial_structure.trial_count)
        warning('trial count and trial info lines mismatch, no trial info attached to output');
    else
        trial_structure.fields = trial_info.fields;
        for i = 1 : trial_info.trial_count
            trial_structure.trials(i).info = trial_info.content(i, :);
        end
    end
end
return;

function ok = has_code(event_group, code)
    ok = false;
    for i = 1 : numel(event_group)
        if (event_group(i) == code)
            ok = true;
        end
    end
return;
