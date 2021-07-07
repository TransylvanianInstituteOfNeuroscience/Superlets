function trial_info = epd_load_trial_info(filepath)
%remember '' and ()

trial_info = struct(   'eti_path', 0,      ...
                'trial_count',0,    ...
                'field_count',0,    ...
                'fields',0,         ...
                'content',0);

f = fopen(filepath);
trial_info.eti_path = filepath;

% load trial count
splitline = strsplit(fgetl(f), ',');
trial_info.trial_count = str2double(splitline{2});

% load field count
splitline = strsplit(fgetl(f), ',');
trial_info.field_count = str2double(splitline{2});

% skip empty line
fgetl(f);

trial_info.fields = strsplit(fgetl(f), ',');

% initialize content cell with desired trial and field counts
trial_info.content = cell(trial_info.trial_count, trial_info.field_count);

for i = 1 : trial_info.trial_count
    
    % read line by line and split into fields
    splitline = strsplit(fgetl(f), ',');
    
    % save each field in the content cell
    for j = 1 : trial_info.field_count
        trial_info.content{i, j} = splitline{j}; 
    end
end

fclose(f);

return;


