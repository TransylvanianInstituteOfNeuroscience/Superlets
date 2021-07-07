function event_groups = epd_parse_event_pattern(event_pattern_string)

group_tokens    = split(event_pattern_string, ',');
event_groups    = { };

for i = 1 : numel(group_tokens)
   
    tokens  = split(group_tokens{i}, {'/', '|'});
    codes   = [];
    
    for j = 1 : numel(tokens)
        
        x = str2double(tokens{j});
        if ~isnan(x)
            codes = [codes, x];
        end
    end
    
    if ~isempty(codes)
        event_groups = [event_groups, codes];
    end
end
        
        