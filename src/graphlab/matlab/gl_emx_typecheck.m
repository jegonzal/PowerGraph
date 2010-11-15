% perform a typecheck of d to ensure that EMLC can translate d
% and that my mxArray -> emxArray converters can understand datatypes.
% This function also changes all arrays to emlcoder.egs dyanamic arrays
% Returns status = 0 on failure, and status = 1 on success
%
% gen returns code when can be used to instantiate a variable of this type
% in an emx m file
function [d, status, gen] = gl_emx_typecheck(d, genprefix)
    if (~exist('genprefix','var'))
        genprefix = 'v';
    end
    gen = '';
    sid = genprefix;
    cname = class(d);
    switch(cname)
        
        
        case {'numeric','integer','int8','uint8','int16','uint16','int32','uint32','int64','uint64','float','single','double'}
            % standard numeric
            % we support. Now, is this a scalar or a matrix.
            if (isscalar(d))
                % scalar. we are done here
                status = 1;
                if (strcmp(cname, 'numeric'))
                    gen = [sid ' = 0.0;'];
                elseif (strcmp(cname, 'integer'))
                    gen = [sid ' = int64(0);'];
                else
                    gen = [sid ' = ' cname '(0);'];
                end
                return;
            else 
                % vector/matrix. convert to dynamic
                d = emlcoder.egs(d(1), size(d) * Inf);
                status = 1;
                if (strcmp(cname, 'numeric'))
                    gen = [sid ' = [0.0];'];
                elseif (strcmp(cname, 'integer'))
                    gen = [sid ' = [int64(0)];'];
                else
                    gen = [sid ' = [' cname '(0)];'];
                end
                gen = [gen '\n' 'eml.varsize(''' sid ''');'];
                return;
            end
            
            
        case {'char'}
            % we only support standard char arrays. not char matrices
            if (isvector(d))
                d = emlcoder.egs('a', [1, Inf]);
                gen = [sid ' = '''';'];
                gen = [gen '\n' 'eml.varsize(''' sid ''');'];
                status = 1;
                return;
            else
                status = 0;
                return;
            end
            
            
        case {'struct'}
            % now structs are complicated. this could be either one struct
            % or a struct array. pull out one element of the struct and
            % parse it recursively
            
            % if it is only one struct, I can fill it in directly
            % otherwise, I will need to go through an intermediary struct
            recursename = '';
            if (numel(d) == 1)
                recursename = sid;
            else
                recursename = [sid '_struct'];
                recursename(find(recursename == '.')) = '_';
            end
            temp = d(1);
            fnames = fieldnames(temp);
            for i = 1:length(fnames)
                % recursively check the field
                [temp.(fnames{i}), status, gen2] = gl_emx_typecheck(temp.(fnames{i}),[recursename '.' fnames{i}]);
                % if fail, return immediately
                gen = [gen '\n' gen2];
                if (status == 0) 
                    return
                end
            end
            % is this one struct or a struct array?
            if (numel(d) == 1)
                % one struct
                d = temp;
            else 
                % vector/matrix!
                d = emlcoder.egs(temp, size(d) * Inf);
                gen = [gen '\n' sid ' = [' recursename '];'];
                gen = [gen '\n' 'eml.varsize(''' sid ''');'];
            end
            status = 1;
        case {'function_handle','logical','cell'}
            % not supported!
            status = 0;
            return;
        otherwise
            % not supported!
            status = 0;
            return;
    end
    
end