function [ res, res_crit_weight, res_values ] = smaowa( input_matrix , delta, wik_ant)
%SMAOWA Implementation of Selective Majority Additive Ordered Weighted 
%   Aritmetic (SMA-OWA) method.
%
%   INPUTS:
%    input_matrix: Matrix where input_matrix(1,:) indicates values and
%       input_matrix(2,:) indicates cardinalities for these values.
%    delta: Cardinality relevance factor used for SMA-OWA.
%    wik_ant (optional): Used for recursive calls. Vector containing the
%       weights calculated in the last step.
%
%   OUTPUTS:
%   res: Final value obtained after applying the SMA-OWA operator to de
%       'input_vector' with a cardinality relevance factor 'delta'.
%   res_crit_weight: Vector containing the weights obtained for each 
%       criterion using the SMA-OWA operator.
%   res_value: Vector where res_value(i) is the value whose weight is
%       res_crit_weight(i).
%    
%
%

    if (nargin < 3)        
        wik_ant = ones(1,size(input_matrix,2))/size(input_matrix,2);
    end   
    
    [values, remainder] = decrease_each_value(input_matrix);
    
    yik = deltavector(delta,values,input_matrix(1,:));
    
    uk = 1 + sum(yik);
    
    wik = (yik + wik_ant)/uk; %calculate w(i,k) using w(i,k-1)
    
    if(not(isempty(values)))
        %RECURSIVE CALL:
      [res, res_crit_weight, res_values] = smaowa(remainder, delta, wik);
    else
        %If there are no values left, the weights are the last calculated ones.
      res = sum(input_matrix(1,:).*wik_ant);
      res_crit_weight = wik_ant;
      res_values = input_matrix(1,:);
    end
    
end

function [values, decreased_matrix] = decrease_each_value(input_matrix)
    % Decreases by 1 each value with cardinality greater than 0
    %
    % INPUT: 
    % input_matrix: Matrix where input_matrix(1,:) are the values and
    %   input_matrix(2,:) are their corresponding cardinalities.
    %
    % OUTPUT:
    % values: Values whose cardinality in the input matrix were greater
    %   than 0 and were decreased by this method.
    % decreased_matrix: input matrix with cardinalities greater than 0 were
    %   decreased by 1.
    %
    values = [];
    decreased_matrix = input_matrix;    
    for i = 1:size(input_matrix,2)
        if(input_matrix(2,i) <= 0)
            decreased_matrix(2,i) = 0;
        else
            decreased_matrix(2,i) = decreased_matrix(2,i) - 1;
            values = [values input_matrix(1,i)];
        end
    end    
end



function [values, remainder] = extract_unique_values(vector)
    % Extracts one instance of each unique value from the vector.
    %
    % INPUTS:
    % vector: vector from where the values are to be extracted.
    %
    % 
    %
    % OUTPUTS: 
    % values: vector containing each value from the input vector but
    % without repetitions.
    % 
    % remainder: input vector taking out one of each different values.    %
    %
    [values, temp_values_pos] = unique(vector);
    for i = temp_values_pos
        vector(i) = [];
    end
    remainder = vector;
end

function [res] = deltavector(delta,values,all_values)
    temp1 = delta*ismember(all_values,values);
    temp2 = (1-delta)*not(ismember(all_values,values));
    res = temp1 + temp2;
end

function [res] = not_in(a,b)
    % Returns unique elements in a that are not in b
    res = setxor(a,b(ismember(b,a)));
end

