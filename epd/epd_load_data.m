function data = epd_load_data(epd, window, markers, channel_indices)

if (~numel(window) == 2)
    error('the window should be specified as two time periods (in seconds) before and after an event marker');
end

if (~isfield(epd, 'version'))
    error('the input structure is not an EPD structure');
end

if (isempty(markers))
    error('marker array empty');
end

if (isstring(channel_indices))
    channel_indices = epd_parse_channel_string(epd, channel_indices);
end
if (isempty(channel_indices))
    error('channel array empty');
end

% 
frame_size      = fix(sum(window) * epd.sampling_rate);
buffer_count    = numel(markers) * numel(channel_indices);
data            = zeros(buffer_count, frame_size);
read_counter    = 1;

frame_before    = fix(window(1) * epd.sampling_rate);
frame_after     = fix(window(2) * epd.sampling_rate);
total_samples   = epd.sample_count;


% 
for i_ch = 1 : numel(channel_indices)
    
    f_ch = fopen([epd.directory filesep epd.channel_files{channel_indices(i_ch)}]);
    if (f_ch == 0)
        error('could not find a channel file');
    end
    
    for i_marker = 1 : numel(markers)
        
        read_begin  = markers(i_marker).timestamp - frame_before;
        read_end    = markers(i_marker).timestamp - frame_after;

        if (read_begin >= 0 && read_end < total_samples)
            fseek(f_ch, read_begin * 4, -1);
            data(read_counter, :) = fread(f_ch, frame_size, 'float');
        else
            % handle edge cases
            if (read_begin < 0 && read_end < total_samples)
                
                % handle begin cutoff
                read_begin = abs(read_begin);
                if (read_begin > frame_size)
                    read_counter = read_counter + 1;
                    continue;
                end
                
                fseek(f_ch, 0, -1);
                data(read_counter, read_begin + 1 : end) = fread(f_ch, read_end - read_begin, 'float');
                
            else
                if (read_begin > 0 && read_end > total_samples)
                    read_end = read_end - total_samples;
                    if (read_end > read_length)
                        read_counter = read_counter + 1;
                        continue;
                    end
                    
                    fseek(f_ch, read_begin, -1);
                    data(read_counter, 1 : frame_size - read_end) = fread(f_ch, frame_size - read_end, 'float');
                else
                    read_begin  = abs(read_begin);  
                    read_end    = read_end - total_samples;

                    fseek(f_ch, read_begin, -1);
                    data(read_counter, read_begin : frame_size - read_end) = ...
                        fread(f_ch, frame_size - read_begin - read_end, 'float');
                end
            end
        end
        
        read_counter = read_counter + 1;
    end
end

return;







