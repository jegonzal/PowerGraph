function generate_update_function_list(updatefunctions)
    % define the macro containing all the update function names
    s = ['#define __UPDATE_FUNCTIONS__ (' num2str(length(updatefunctions)) ',('];
    for i = 1:length(updatefunctions)
        s = [s updatefunctions{i}];
        if (i ~= length(updatefunctions))
            s = [s ','];
        end
    end
    s = [s '))\n'];
    fprintf(s);
    
    
    f = fopen('update_function_array.hpp', 'w');
    % standard HPP headers
    fprintf(f, '#ifndef UPDATE_FUNCTION_ARRAY_HPP\n');
    fprintf(f, '#define UPDATE_FUNCTION_ARRAY_HPP\n');
    % include all the update function headers
    for i = 1:length(updatefunctions)
        fprintf(f, ['#include "' updatefunctions{i} '.h"\n']);
    end
    % output the UPDATE_FUNCTIONS macro
    fprintf(f, s);
    fprintf(f, '#endif\n');
    fclose(f);
end