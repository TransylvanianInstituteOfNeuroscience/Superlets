
dataset_folder  = 'c:\_data\eeg\dots_30\dots_30_001';
dataset_name    = 'dots_30_001';

% load the dataset and the trial info
epd = epd_open([dataset_folder filesep dataset_name '.epd']);
eti = epd_load_trial_info([dataset_folder filesep dataset_name '.eti']);

% create a trial structure
trials = epd_parse_trial_structure(epd, '128,150,129,1/2/3,131,132', eti);
trials_f = epd_filter_trials(trials, 'ResponseID', '3');

% get markers and load data
markers = epd_get_markers(trials, '1/2/3');

% get data
data = epd_load_data(epd, [1, 1], markers, [ 2 ]);