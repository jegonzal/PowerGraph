%% Construct a discrete table factor
%   
%   factor = table_factor(vars, logP)
%
% vars: array of variable ids (e.g., [1,2,4] )
% logP: tensor representing the log potential values (e.g., ones(3,7,2)
%    where variable 1 takes on 3 states variable 2 takes on 7 states and
%    variable 4 takes on 2 states.
%   
% A table factor represents a factor or potential over a small set of
% discrete variables.  Forexample if we wanted to encode a similarity
% funciton over a pair of variables x1 and x2 we could define the table
% factor:
%
%    psi(x1,x2) = exp( |x1 - x2| )
% 
% Assuming x1 and x2 take on 4 and 3 values respectively we could define
% the matrix (table) representation of psi:
%    
%              0  1  2 
%  tbl = exp(  1  0  1  )
%              2  1  0
%              3  2  1
% 
% We can build a table factor representing this as:
%
%  factor = table_factor([1, 2], tbl);
%
function factor = table_factor(vars, data) 
factor.vars = sort(uint32(vars));
factor.logP = double(data);
end

%%
% Originally I had hoped to used a matlab class but unfortunatley mex
% support for classes is limited resulting in substantial performance
% penalties when accessing fields.

% classdef table_factor
%   properties (SetAccess = private)
%     variables;
%   end
%   properties
%     logP;
%   end
%   methods
%     %% the variables should be a d dimensional array
%     function obj = table_factor(vars, data)
%       obj.variables = uint32(vars);
%       obj.logP = double(data);
%     end
%     
%   end
%   
% end

