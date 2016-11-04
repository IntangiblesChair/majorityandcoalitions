function [ res ] = choquet_c( values, perm_map )
    % Calculates the choquet integral for the values given using the 
    % weights stored in a permutation map.
    %
    % INPUTS: 
    % values: vector of values obtained evaluating a given alternative in
    % each criteria.
    % perm_map: map containing the combinations of criteria an their respective 
    % weigths.
    %
    % OUTPUTS:
    % res: Choquet integral for 'values' using the weigths given in
    %  'perm_map'
    %

    res = 0;
    keys = char(65):char(64+length(values));

    while max(values)>0
        pos = find(values == min(values));
        pos = pos(1);
        this_value = values(pos);
        res = res + perm_map(keys)*this_value;
        values = [values(1:pos-1) values(pos+1:end)];
        values = values-this_value;
        keys = [keys(1:pos-1) keys(pos+1:end)];        
    end    
end