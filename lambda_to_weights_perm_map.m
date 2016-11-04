function [ res_map ] = lambda_to_weights_perm_map( weights_vector, coalitions_lambda_map )
    % Generates a map with the weigths of the combinations taking into
    % account criteria coalitions.
    %
    % INPUTS:
    % weights_vector: vector containing the weights calculated for each
    % criteria. It is assumed that sum(weigths_vector) == 1
    %
    % coalitions_lambda_map: Map containing the keys for each criteria coalition 
    %  in the form key='Tag' value='Value' where 'Tag' is a string specifying 
    %  which criteria participate in the criterion and 'Value' is the lambda 
    %  calculated for the criteria coalition.
    %  As an example, {'AC','-2.47'} indicates a lambda value of -2.47 for the 
    %  coalition formed between the first and third criterions in the order given 
    %  by the weights_vector.
    %
    % OUTPUTS: 
    % res_map: Map with the weights of the permutations of weight_vector
    %  taking into account criteria coalitions. Normalized between [0,1].
    
    %
    % Since this algorithms takes into account criteria coalitions, each
    % combination has to be calculated once for each permutation leading to
    % it. This is because if [AB] is influenced by a criteria coalition, 
    % [AB]+[C] is not neccesarily equal to [A]+[BC].
    %


    keys = char(65):char(64+length(weights_vector));
    keys = num2cell(keys);
    values = num2cell(weights_vector);
    step_map = containers.Map(keys,values); %Contains the keys calculated in the previous step
    res_map = containers.Map(); %Final solution containing every combination.
    if nargin < 2
        coalitions_lambda_map = containers.Map();
    end
    
    
    for i = 1:length(weights_vector)
        temp_map = containers.Map(); %Map calculated in this iteration of i. Will contain every combination with length(i+1).
        for j = step_map.keys
            for k = 1:length(weights_vector)
                this_char = char(k+64);
                if isempty(strfind(char(j),char(k+64)))
                    this_key = sort([char(j) this_char]); %Key composition, ie: 'AC' + 'B' = 'ABC'
                    
                    %Obtaining lambda
                    if nargin > 1
                        if isKey(coalitions_lambda_map, this_key) %If the coalition exists in coalitions_lambda_map, calculate lambda.
                            lambda = coalitions_lambda_map(this_key);
                        else
                            lambda = 0; %If there is no coalition, lambda = 0.
                        end
                    else
                        lambda = 0; %If no coalitions_lambda_map is given, it is assumed that there are no coalitions.
                    end 
                    
                    %Weight calculation.
                    this_value = step_map(char(j)) + weights_vector(k) + lambda*step_map(char(j))*weights_vector(k);
                    
                    %Weight storage in the map.
                    if isKey(temp_map,this_key)
                        %If a value already exists, only update if the new
                        %value is greater than the old value.
                        if this_value > temp_map(this_key)
                            temp_map(this_key) = this_value;
                        end
                    else
                        %If there was no value stored, it gets stored now.
                        new_map = containers.Map(num2cell(this_key,2),num2cell(this_value));                        
                        temp_map = [temp_map; new_map];
                    end                    
                end
            end %for k
        end %for j
        res_map = [res_map; step_map];
        step_map = temp_map;
    end %for i

    res_map = normalize(res_map);

end

function [map] = normalize (map)
    max_val = max(cell2mat(map.values));
    for i = map.keys
        map(char(i)) = map(char(i))./max_val;
    end
end

