function [ res, perm_map] = coalition_and_mayority( experts_criteria_labels, alternatives_file, criteria_file, coalitions, coalition_linguistic_values)
    % Implements Selective Majority Additive OWA with Criteria Coalitions
    %
    % INPUTS: 
    % experts_criteria_labels: Vector containing the values of each linguistic label.
    %
    % alternatives_file: Path of a file containing values obtained 
    %  for every alternative in each criteria.
    %
    % criteria_file: Path of a file containing a criteria matrix where
    %  element (i,j) represents the cardinality of criteria i for value j.
    %
    % coalitions: Path of a file containing coalition data or map containing 
    %  the keys for each criteria coalition  in the form key='Tag' 
    %  value='Value' where tag is a string specifying which criteria 
    %  participate in the criterion and value is a number in 
    %  the interval [-1,1] indicating it's sign and weight.
    %
    % coalition_linguistic_values: Values of the linguistic tags for
    % criteria coalitions.
    %
    % OUTPUTS: 
    % res: Matrix where res(i,j) is the calculated value for the
    %  criteria j using SMA-OWA and criteria coalition.
    %
    % perm_map: containers.Map with the weights of every single
    %  combination of criteria.
    %
    % NEEDS:
    % lambda_to_weights_perm_map: Calculate perm_map.
    % smaowa: Implementation of the SMA-OWA Algorithm.
    % choquet_c: Implementation of Choquet integral using perm_map.
    
    
    %========= GLOBAL SETTINGS ==============
    ZERO_TOLERANCE = 1e-10; %Values smaller than this will be assumed to be zero.  
    
    
    %========= READING INPUT DATA ============
    %Read the experts criteria file
    num_label_values = length(experts_criteria_labels);
    if ischar(criteria_file)
        fi = fopen(criteria_file,'r'); 
        temp = fscanf(fi,'%f',[double(num_label_values) inf]);     
        fclose(fi);    
        cardinalities = temp';
    else
        %If not a file, use the given matrix.
        cardinalities = criteria_file;
    end
    num_criteria = size(cardinalities,1);
    
    
    %Read the alternative evaluations file
    if ischar(alternatives_file)
        fi = fopen(alternatives_file,'r'); 
        temp = fscanf(fi,'%f',[double(num_criteria) inf]); 
        fclose(fi);
        A = temp';
    else
        %If not a file, use the given matrix.
        A = alternatives_file;
    end
    
    %Read the coalition file, when available.
    if nargin < 5
        %When no coalitions are given, assume that there are none.
        coalition_linguistic_values = 0;
        coal_matrix = zeros((num_criteria+1)*(num_criteria/2)-num_criteria,length(coalition_linguistic_values));
    else
        if ischar(coalitions)
            %When a filepath is given generate coalition map from file
            fi = fopen(coalitions,'r');
            num_posib_lv = length(coalition_linguistic_values);
            temp = fscanf(fi,'%f',[double(num_posib_lv) inf]);
            fclose(fi);    
            coal_matrix = temp';
        else
            %If a coalition matrix is already provided, use it.
            coal_matrix = coalitions;
        end
    end
    
    
    
    %========= PROCESSING DATA ============
    
    %---- STEP 1: Use SMA-OWA to calculate weights -----
    for i=1:num_criteria
        delta = 1-1/(2+std(cardinalities(i,:))); %Use cardinalities's dispersion to calculate delta.
        smaowa_value = smaowa([experts_criteria_labels; cardinalities(i,:)],delta);
        W(i,:) = smaowa_value;
    end
    W = W';
    W = W./sum(W); %Normalize weights
    
    
    
    %---- STEP2: Calculate coalitions -----
    coalition_lambda_map = containers.Map; %At the end of Step 2 will contain the lambda value for each coalition.
    iterator = 1;        
    for i = 1:length(experts_criteria_labels)
        for j = (i+1):length(experts_criteria_labels)
            
            this_iteration_cardinalities = coal_matrix(iterator,:);
            
            %Pick every value with a cardinality greater than zero.
            non_zero_linguistic_values = [coalition_linguistic_values(find(abs(this_iteration_cardinalities)>ZERO_TOLERANCE)); this_iteration_cardinalities(this_iteration_cardinalities ~= 0)];
            
            %Calculate lambda.
            if isempty(non_zero_linguistic_values)
                %If no criteria coalition is present, lambda = 0.
                lambda = 0;
            else
                %Calculate lambda_max and lambda_range
                if(abs(W(i)*W(j)) > ZERO_TOLERANCE)
                    lambda_min = (max(W(i),W(j))-(W(i)+W(j)))/(W(i)*W(j));
                else
                    lambda_min = 0;
                end
                lambda_max = abs(lambda_min);
                lambda_range = lambda_max/((length(coalition_linguistic_values)-1)/2);
                delta = 1-1/(2+std(coal_matrix(iterator,:)));  %Use cardinalities' dispersion to calculate delta.
                
                
                %Use SMA-OWA to calculate lambda
                smaowa(non_zero_linguistic_values,delta);
                lambda = smaowa(non_zero_linguistic_values,delta)*lambda_range;
            end
            
           
            %If lambda ~= 0, store the lambda value for this coalition.
            if(abs(lambda) > ZERO_TOLERANCE)
                temp_map = containers.Map(num2cell([num2char(i) num2char(j)],2),lambda);
                coalition_lambda_map = [coalition_lambda_map ; temp_map];
            end
            
            iterator = iterator +1;
        end
    end    
    
    %---- STEP 3: Agregate with Choquet each alternative -----
    perm_map = lambda_to_weights_perm_map(W,coalition_lambda_map); %Obtain the permutation map.
    
    res = [];
    for i=1:size(A,1)        
        res_alt = choquet_c(A(i,:),perm_map);  %Use choquet integral to aggregate
        res = [res res_alt];
    end

end



function [res] = num2char(num)
    %Transforms a number to char 1->A, 2->B and so on.    
    res  = char(num+64);
end
    
    
    
    
    
    
    
    
    
    
    
    


