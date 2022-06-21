function epd = epd_open(path)

epd = struct('version',                 0,  ... 
            'directory',                0,  ...
            'channel_count',            0,  ...
            'sampling_rate',            0,  ...
            'sample_count',             0,  ... 
            'channel_files',            0,  ...
            'event_timestamps_file',     0,  ...
            'event_codes_file',         0,  ... 
            'event_count',              0,  ...
            'channel_names',            0,  ... 
            'average_reference_count',  0,  ... 
            'averagerefs',              0);
f = fopen(path);

[epd.directory] = fileparts(path);

line = fgetl(f);
while ischar(line)
    switch line
        case 'Format version:'
            line = fgetl(f);
            epd.version = str2double(line);
            
        case 'Number of EEG channels:'
            line = fgetl(f);
            epd.channel_count = str2double(line);
            
        case 'Sampling frequency (Hz):'
            line = fgetl(f);
            epd.sampling_rate = str2double(line);
            
        case 'Total number of samples:'
            line = fgetl(f);
            epd.sample_count = str2double(line);
            
        case 'List with filenames that hold individual channel samples (32 bit IEEE 754-1985, single precision floating point; amplitudes are measured in uV):'
            epd.channel_files = cell(epd.channel_count, 1);
            for i = 1:epd.channel_count
                line = fgetl(f);
                epd.channel_files{i} = line;
            end
            
        case 'File holding event timestamps; timestamp is in samples; (32 bit signed integer file):'
            line = fgetl(f);
            epd.event_timestamps_file = line;
            
        case 'File holding codes of events corresponding to the event timestamps file; timestamp is in samples; (32 bit signed integer file):'
            line = fgetl(f);
            epd.event_codes_file = line;
            
        case 'Number of events (size of the list with event timestamps and list with event codes):'
            line = fgetl(f);
            epd.event_count = str2double(line);
            
        case 'List with labels of EEG channels:'
            epd.channel_names = cell(epd.channel_count, 1);
            for i = 1:epd.channel_count
                line = fgetl(f);
                epd.channel_names{i} = line;
            end
            
        case 'Number of channels that were considered for computing an average reference:'
            line = fgetl(f);
            epd.average_reference_count = str2double(line);
            
        case 'List of channels (index name) that were taken into computing an average reference (if empty, then signals were not referenced when importing to EPD; be careful with Biosemi datasets that are recorded unreferenced...):'
            if epd.average_reference_count ~= 0
                epd.averagerefs = cell(epd.average_reference_count, 1);
                for i = 1:epd.average_reference_count
                    line = fgetl(f);
                    epd.averagerefs{i} = line;
                end
            end
    end
        
    line = fgetl(f);
    
end

fclose(f);

   


end