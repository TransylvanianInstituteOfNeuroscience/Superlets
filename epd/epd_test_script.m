
dataset_folder  = 'c:\_data\eeg\dots_30\dots_30_001';
dataset_name    = 'dots_30_001';
trial_structure = '128,150,129,1/2/3,131,132';

% load the dataset and the trial info
epd = epd_open([dataset_folder filesep dataset_name '.epd']);
eti = epd_load_trial_info([dataset_folder filesep dataset_name '.eti']);

% create a trial structure
trials = epd_parse_trial_structure(epd, trial_structure, eti);
trials_f = epd_filter_trials(trials, 'ResponseID', '3');

% get markers and load data
markers = epd_get_markers(trials, '1/2/3');

% get data
seconds_before  = 1;
seconds_after   = 2;
channel_indices = [ 2 ];

data = epd_load_data(epd, [seconds_before, seconds_after], markers, channel_indices);