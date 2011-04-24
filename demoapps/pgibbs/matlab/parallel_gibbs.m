%% This file represents a wrapper to the parallel Gibbs sampling tools

function [samples, marginals] = ...
   parallel_gibbs(factors, nsamples, frequency, method)

if(~strcmp(class(factors), 'cell'))
   error('The factors argument must be a cell array of table_factors');
end


max_var = 0;
% Check the factors data structure
for i = 1:length(factors)
   if(~strcmp(class(factors{i}), 'table_factor'))
      disp(factors{i});
      error('Factor %d is not of type table_factor', i);
   end
   % Get the maximum variables
   max_var = max(max(factors{i}.variables(:)), max_var);
   if(min(factors{i}.variables) <= 0) 
      disp(factors{i});
      error('Factor %d has 0 valued variables', i);
   end   
end

% check all the variables have consistent sizes;
vars = 1:max_var;
var_sizes = zeros(max_var, 1);
for i = 1:length(factors)
   current_sizes = var_sizes(factors{i}.variables);
   dims = size(factors{i}.logP)';
   dims = dims(dims > 1);
   if(length(dims) ~= length(current_sizes))
      error(['The number of variables %d in factor %d does not match ' ...
         'the number of dimensions %d.'], ...
         length(current_sizes), i, length(dims));
   end   
   ind = current_sizes > 0;
   % all elements in ind have been set and should just match
   if( ~isempty(find(current_sizes(ind) ~= dims(ind), 1)) ) 
      errorind = find(current_sizes(ind) ~= dims(ind), 1);
      error(['Variable %d has already been seen having size %d ', ...
         'but was just now observed to have size %d in factor %d.'], ...
          factors{i}.variables(ind(errorind)), ...
          current_sizes(ind(errorind)), ...
          dims(ind(errorind)), ...
          i);
   end
   var_sizes(factors{i}.variables(~ind)) = dims(~ind);
end
unset_vars = find(var_sizes(:) == 0);
if(~isempty(unset_vars)) 
   error(['The following variables were not set correctly: ', ...
      mat2str(unset_vars)]);   
end


chromatic_gibbs(factors);



end