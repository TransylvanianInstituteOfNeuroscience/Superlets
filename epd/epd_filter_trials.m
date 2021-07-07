function new_trials = epd_filter_trials(trial_structure, varargin)

if (isempty(trial_structure.fields))
    error('trial structure contains no trial info');
end

narginchk(2, inf);
criteria_count = fix(numel(varargin) / 2);
if (rem(numel(varargin), 2) ~= 0)
    warning('possible error detected in varargin: uneven number of keys and values specified');
end


new_trials              = [];
new_trials.trials       = [];
new_trials.trial_count  = 0;
new_trials.fields       = trial_structure.fields;

for i = 1 : trial_structure.trial_count
    
    keep = true;
    
    for i_arg = 1 : criteria_count
        field_idx = find(strcmp(trial_structure.fields, varargin{(i_arg - 1) * 2 + 1}));
        if (~isempty(field_idx) && ~strcmp(trial_structure.trials(i).info{field_idx}, varargin{i_arg * 2}))
            keep = false;
            break;
        end
    end
    
    if (keep)
        new_trials.trials = [new_trials.trials, trial_structure.trials(i)];
    end       
end

new_trials.trial_count = numel(new_trials.trials);

return;


