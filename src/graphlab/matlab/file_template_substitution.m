% function file_template_substitution(TEMPLATE, TARGET, SUBST) 
%
% Copies TEMPLATE to TARGET, substituting a few macros as described by
% SUBST
%
%   TEMPLATE: File name for the source template file
%   TARGET: File name for the source template file
%   SUBST: A struct where each field is a string.
%          The field name ->string mapping is the substitution rule
%
%  For each field FIELD with value VALUE in the struct STRUCT.
%  All occurances of '%#FIELD#%' within the TEMPLATE file will be replaced
%  with the VALUE
function file_template_substitution(templatefile, targetfile, substs)
    srcfile = fopen(templatefile);
    destfile= fopen(targetfile, 'w');
    
    fnames = fieldnames(substs);
    
    while(1)
        line = fgets(srcfile);
        if (line == -1) 
            break
        end
        for i = 1:length(fnames)
            line = strrep(line, ['%#' fnames{i} '#%'], substs.(fnames{i}));
        end
        fprintf(destfile,'%s',line);
    end
    fclose(destfile);
    fclose(srcfile);
end