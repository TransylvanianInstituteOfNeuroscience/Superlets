function markers = epd_load_events(epd)

markers = [];

f_evc = fopen([epd.directory filesep epd.event_codes_file]);
f_evt = fopen([epd.directory filesep epd.event_timestamps_file]);

if (f_evc == 0 || f_evt == 0)
    error('could not open event codes or timestamps file');
end

evc = fread(f_evc, inf, 'int32');
evt = fread(f_evt, inf, 'int32');
fclose(f_evc);
fclose(f_evt);

if (numel(evc) < epd.event_count || numel(evt) < epd.event_count)
    error('number of read events is smaller than the number specified in the epd');
end

for i = 1 : epd.event_count
    markers(i).code         = evc(i);
    markers(i).timestamp    = evt(i);
end

return;
