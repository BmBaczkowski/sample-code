function T = get_prolific_data(prolificfile, keyfile)

prolific    = readtable(prolificfile, 'PreserveVariableNames', 1);
txt         = fileread(keyfile);     
key         = jsondecode(txt);

data        = prolific(:, {'participant_id', 'age', 'Sex'});
pids        = fieldnames(key);

sz          = [length(pids), 3];
varTypes    = {'double','double','string'};
T           = table('Size',sz,'VariableTypes',varTypes, 'VariableNames', {'subject', 'age', 'sex'});

for i = 1:length(pids)
   
    pid     = pids{i};
    indx    = find(ismember(data{:,1},pid(2:end)));
    T{i,1}  = str2num(key.(pid));
    T{i,2}  = data{indx,'age'};
    T{i,3}  = data{indx,'Sex'};
    
end

end