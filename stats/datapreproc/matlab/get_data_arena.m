function [data] = get_data_arena(sourcedir, prolificdata, qcdata)

data = [];
for subject = prolificdata.subject'
  
    fprintf('Subject: %03d\n', subject);
    sourcefile_coords    = fullfile(sprintf('%s/%03d/arena_coords.csv', sourcedir, subject));
    sourcefile_graph     = fullfile(sprintf('%s/%03d/arena.csv', sourcedir, subject));
    
    if exist(sourcefile_graph, 'file')
        coords   = readtable(sourcefile_coords);
        graph    = readtable(sourcefile_graph);

        % normalize
        graph.x = (graph.x - coords.left) / (coords.width);
        graph.y = (graph.y - coords.top) / (coords.height);
        
        sz = [size(graph,1) 3];
        T2 = table('Size', sz, 'VariableTypes', {'double','double', 'string'}, 'VariableNames', {'subject', 'age', 'sex'});
        indx = prolificdata.subject == subject;
        T2{:,1} = repmat(prolificdata.subject(indx), size(graph,1), 1);
        T2{:,2} = repmat(prolificdata.age(indx), size(graph,1), 1);
        T2{:,3} = repmat(prolificdata.sex(indx), size(graph,1), 1);
        
        T = [T2 graph];
        data = [data; T];
    end
end

end