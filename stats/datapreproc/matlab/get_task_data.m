function data = get_task_data(sourcedir, prolificdata, taskduration)

data = [];
for subject = prolificdata.subject'
    
    
    fprintf('Subject: %03d\n', subject);
    for block = 1:11
        sourcefile = fullfile( sprintf('%s/%03d/block%02d.csv', sourcedir, subject, block));
        if exist(sourcefile, 'file')
            
            T1      = readtable(sourcefile);
            sz      = [size(T1,1) 5];
            T2      = table('Size', sz, 'VariableTypes', {'double','double', 'string', 'double', 'double'}, 'VariableNames', {'subject', 'age', 'sex', 'block', 'duration'});
            indx    = prolificdata.subject == subject;
            T2{:,1} = repmat(prolificdata.subject(indx), size(T1,1), 1);
            T2{:,2} = repmat(prolificdata.age(indx), size(T1,1), 1);
            T2{:,3} = repmat(prolificdata.sex(indx), size(T1,1), 1);
            T2{:,4} = repmat(block, size(T1,1), 1);
            indx    = taskduration(:,1) == subject & taskduration(:,2) == block;
            T2{:,5} = repmat(taskduration(indx,3), size(T1,1), 1);
            T       = [T2 T1];
            data    = [data; T];

        end
    end

end