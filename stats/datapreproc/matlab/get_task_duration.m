function data = get_task_duration(jatosfile, keyfile)

txt         = fileread(keyfile);     
key         = jsondecode(txt);
pids        = fieldnames(key);

jatos       = readtable(jatosfile, 'Delimiter','comma');

data = [];
for i = 1:length(pids)
    
    pid     = pids{i};
    sub     = key.(pid);
    sub     = str2num(sub);
    pid     = pid(2:end);
    
    offsets = [strcmp('Offset', jatos{:,1}), ...
                strcmp('task', jatos{:,2}), ...
                strcmp('seq', jatos{:,3}), ...
                strcmp(pid, jatos{:,5})];
    offsets = all(offsets, 2);
    offsets = find(offsets);
    
    block = 0;
    for offset = offsets'
        block = block + 1;
        off = jatos{offset,6};
        onsets = strcmp('Onset', jatos{1:offset,1}) & ...
                strcmp('task', jatos{1:offset,2}) & ...
                strcmp('seq', jatos{1:offset,3}) & ...
                strcmp(pid, jatos{1:offset,5});
        onsets = find(onsets);
        onset = onsets(end);
        on  = jatos{onset,6};
        duration = off - on;
        
        data = [data; [sub block duration]];
    end
    
   
end
end