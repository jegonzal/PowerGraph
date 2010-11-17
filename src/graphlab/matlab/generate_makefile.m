% generate a Makefile which builds a mex library. 
% This will detect the locations of the mex headers and mex libraries
% automatically. The g++ compiler will be used. Since this uses mex 
% configuration information to generate the Makefile, mex setup should
% be done first, selected gcc/g++ as the compiler.
%
% makefilename : The make file is created with this filename
% mexbasename : The name of the mex function. Do not add the extension.
%               This will be autodetected
% cppfiles : All the C/CPP files to compile
% includepaths : directories to search for header files. 
% libpaths : directories to search for lib files
% cppflags : Additional flags to pass to the CPP compiler
% ldflags : additional flags to pass to the linker
function generate_makefile(makefilename, binname, cppfiles, includepaths, libpaths, cppflags, ldflags, ismex, objprefix)
    if (nargin < 8)
        disp('insufficient arguments');
        return;
    end
    if (~exist('objprefix', 'var')) 
        objprefix = '';
    end

    % get my current path
    [st,~] = dbstack('-completenames');
    mfiledirectory = st.file;
    % get the graphlab/matlab directory
    slashes = strfind(mfiledirectory, '/');
    glmatlabdir = mfiledirectory(1:(slashes(end)-1));
    
    % get the makefile location
    slashes = strfind(makefilename, '/');
    makefilepath = '.';
    if (~isempty(slashes))
        makefilepath = makefilename(1:(slashes(end)-1));
    end
    extraopts = {};
    if (ismex)
        binname = [binname '.' mexext];
        % get the matlab compilation options
        disp([glmatlabdir '/get_mex_params.sh ' makefilepath '/config.mk']);
        unix([glmatlabdir '/get_mex_params.sh ' makefilepath '/config.mk']);

        f = fopen([makefilepath '/config.mk']);
        i = 1;
        while (~feof(f))
            extraopts{i} = fgets(f);
            extraopts{i} = strrep(extraopts{i}, '-ansi', '');
            i = i + 1;
        end
        fclose(f);

    end
    % add the include paths into the cppflags
    for i = 1:length(includepaths)
        cppflags = ['-I"' includepaths{i} '" ' cppflags];
    end
    % add the libpaths into the ldflags
    for i = 1:length(libpaths)
        ldflags = [ '-L"' libpaths{i} '" ' ldflags];
    end
    f = fopen(makefilename, 'w');

    fprintf(f, 'default_target: %s\n\n', binname);
    % output the extra options
    for i = 1:length(extraopts)
        fprintf(f, '%s\n', extraopts{i});
    end
    if (ismex) 
        fprintf(f, 'CPPFLAGS = $(MEXCXXFLAGS) %s \n\n', cppflags);
        fprintf(f, 'LDFLAGS = $(MEXLDFLAGS) %s\n\n', ldflags);
    else
        fprintf(f, 'CPPFLAGS = %s \n\n', cppflags);
        fprintf(f, 'LDFLAGS = %s\n\n', ldflags);
    end
    % a target for each object file
    allofiles ='';
    for i = 1:length(cppfiles)
        cppfile = cppfiles{i};
        % turn it into a .o
        dotpos = max(strfind(cppfile, '.'));
        ofile = [objprefix cppfile(1:dotpos), 'o'];
        % make the o file local
        slashpos = strfind(ofile, '/');
        if (~isempty(slashpos))
            ofile = ofile(slashpos(end)+1:end);
        end
        allofiles = [allofiles ' ' ofile];
        
        fprintf(f, '%s: %s Makefile\n', ofile, cppfile);
        fprintf(f, '\tg++ -c $(CPPFLAGS) $< -o $@\n\n\n');
    end
    % final mex target
    fprintf(f, '%s: Makefile %s\n', binname, allofiles);
    fprintf(f, '\tg++ -o $@ %s $(LDFLAGS)\n\n', allofiles);
    
    fprintf(f, 'clean:\n');
    fprintf(f, '\trm -f %s %s', allofiles, binname);
    fclose(f);
    %-Wl,-rpath-link,/usr/local/MATLAB/R2010b/bin/glnxa64 -L/usr/local/MATLAB/R2010b/bin/glnxa64 -lmx -lmex -lmat -lm
end