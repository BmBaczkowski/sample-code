function data = get_distances(arena)

% community matrix of 
% -1 for between
% 1 for within

%% create community matrix
C = ones(11, 11) * -1;


C(1, [2, 3, 4, 5])      = 1;
C(2, [1, 3, 4, 5])      = 1;
C(3, [1, 2, 4, 5])      = 1;
C(4, [1, 2, 3, 5])      = 1;
C(5, [1, 2, 3, 4])      = 1;
C(7, [8, 9, 10, 11])    = 1;
C(8, [7, 9, 10, 11])    = 1;
C(9, [7, 8, 10, 11])    = 1;
C(10, [7, 8, 9, 11])    = 1;
C(11, [7, 8, 9, 10])    = 1;

I = eye(11);
C(logical(I(:))) = 0;
C(6, :) = 0;
C(:, 6) = 0;
C = squareform(C);
%%

data = [];
for subject = unique(arena.subject)'
    
    indx    = arena.subject == subject;
    
    D       = pdist([arena.x(indx) arena.y(indx)]);
    
    d1      = mean(D(C == 1));
    d2      = mean(D(C == -1));
    
    val     = d1 - d2;
    
    sz = [1, 2];
    T  = table('Size', sz, 'VariableTypes', {'double','double'}, 'VariableNames', {'subject', 'dist'});
    T{1, 1} = subject;
    T{1, 2} = val;
    
    data = [data; T];
end



end