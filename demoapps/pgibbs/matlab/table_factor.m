function factor = table_factor(vars, data) 
factor.vars = sort(uint32(vars));
factor.logP = double(data);
end


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

