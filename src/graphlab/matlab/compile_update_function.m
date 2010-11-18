% compile_update_function    Compiles a collection of update functions into
%                            a graphlab program.
%
%   compile_update_function(UPDATE_FNS, EX_VERTEX, EX_EDGE, 
%                           GLLIB_DIR, WORKING_DIR, OUT_NAME,
%                           OPT_LEVEL = 0, PAR_COMPILES = 4, 
%                           NO_TCMALLOC = 0);
%    
%   This function is the wrapper around the Matlab -> GraphLab compilation
%   process. It will use emlc to compile a collection of update function
%   names as listed as strings in UPDATE_FNS into binary executables in
%   the WORKING_DIR.
%
%
%
% Input: 
%
%   UPDATE_FNS: is a cell array of strings, where each string is the name
%               of a GraphLab update function to compile.
%               (see the README for details on the update functions)
% 
%   EX_VERTEX: This should be representative example of the data on 
%              a vertex.
%
%   EX_EDGE  : This should be representative example of the data on 
%              a edge.
%
%   GLLIB_DIR: This should point to a directory where the GraphLab binary
%              libraries can be found. For instance, this could be
%              empty if all graphlab libraries are already in /usr/lib.
%              Or if you compiled GraphLab without installing it, 
%              this could point to the release/src/graphlab directory.
%
%   WORKING_DIR:  All intermediate source files and Makefiles will be 
%                 stored in this directory. The resultant binaries will
%                 also be stored here. This directory will be created
%                 if it does not already exist
%
%   OUT_NAME:  This is the base name to use for the final compiled
%              binaries/m files.
%
%   OPT_LEVEL:  Optimization level. corresponds to the -O[N] flag in 'gcc'.
%               Defaults to 0.
%
%   PAR_COMPILES:  This function will automatically begin building the 
%                  Makefiles. This is the parallelization level for the 
%                  build and corresponds to -j[N] flag in 'make'.
%                  Defaults to 4.
%
%   NO_TCMALLOC:  This function will try to detect the existance of 
%                 libtcmalloc through the use of 'whereis'. If detected,
%                 it will automatically try to link the binaries against
%                 libtcmalloc. If this flag is set to non-zero, we will
%                 never link against tcmalloc even if detected.
%                 Defaults to false.
%
%
%
% Output:
%
%   The script will poroduce a Makefile in the WORKING_DIR which will 
%   compile the following 3 binaries (OUT_NAME)_save_graph,
%   (OUT_NAME)_binary and (OUT_NAME)_load_graph.
%
%
%
%
%
%   (OUT_NAME)_save_graph is a mex library
%
%            (OUT_NAME)_save_graph(VDATA, ADJ, EDATA, ...
%                                       OPTIONS, GRAPHFILE, STRICT)
%
%   This function will serialize the graph described by VDATA, ADJ, EDATA 
%   as well as OPTIONS.initial_schedule to the file GRAPHFILE. If STRICT 
%   is set stricter type checking will be used. 
%
%      VDATA: is a cell array of vertex data where VDATA{i} is the data on 
%             vertex i. The length of this cell array is the number of 
%             vertices in the graph.
%
%      ADJ :  A (sparse) adjacency matrix. If ADJ(i,j) > 0, there is an 
%             edge from vertex i to vertex j, where the data on the edge 
%             i->j is EDATA{ADJ(i,j)} The maximum dimensions of ADJ is 
%             #vertices * #vertices.
%
%      EDATA: A cell array of edge data. The ADJ adjacency matrix 
%             references into this array. It is possible for many edges to 
%             share the same data element. For instance, if I have many 
%             edges with the same value, it is possible for all their ADJ 
%             entries to point to a single entry in EDATA.
%
%      OPTIONS: A struct of only one field "initial_schedule". 
%               (leaving room for future extensions)
%
%         OPTIONS.initial_schedule = [SCHED1, SCHED2 ...]
%                 A struct array describing the initial schedule where each
%                 sched entry has the following format:
%
%         SCHED.update_function : A string. The name of the update function
%         SCHED.vertices        : array of vertex ids to update
%         SCHED.priority        : same size as vertices. The priority for 
%                               each update. Must be >0
%
%         For instance, if I have two update functions 'u1' and 'u2'
%         and I would like u1 to update vertices 1:100 with high priority,
%         and u2 to update vertices 101:200 with lower priority. 
%         I could fill in the OPTIONS:
%
%             OPTIONS.initial_schedule(1).update_function = 'u1'
%             OPTIONS.initial_schedule(1).vertices = 1:100
%             OPTIONS.initial_schedule(1).priority = ones(1,100);
%
%             OPTIONS.initial_schedule(2).update_function = 'u2'
%             OPTIONS.initial_schedule(2).vertices = 101:200
%             OPTIONS.initial_schedule(2).priority = 0.1 * ones(1,100);
%
%      GRAPHFILE: the file to output to.
%
%      STRICT: numeric. Whether strict type checking will be used.
%    
%
%
%
%
%
%   (OUT_NAME)_binary is a standalone program
%
%   This program takes in the standard GraphLab command line options, 
%   (run (OUT_NAME)_binary --help for details), as well two additional
%   options. 
%
%      --ingraphfile=[INGRAPH] --outgraphfile=[OUTGRAPH]
%
%   INGRAPH will read the GRAPHFILE generated by (OUTPUT_M_NAME)_save_graph
%   and run GraphLab on it, using the parameters specified on the command 
%   line, as well as the initial schedule passed into the OPTIONS parameter 
%   of (OUTPUT_M_NAME)_save_graph. When complete the result will be stored
%   in OUTGRAPH which can be read be (OUTPUT_M_NAME)_load_graph.
%
%
%
%
%
%
%   (OUT_NAME)_load_graph is a mex library
%
%            [VDATA, ADJ, EDATA] = (OUT_NAME)_load_graph(GRAPHFILE)
%
%   This function will deserialize the graph in GRAPHFILE and return the 
%   result in VDATA, ADJ, EDATA. See the documentation for 
%   (OUT_NAME)_load_graph for the format of VDATA, ADJ, EDATA.
%  
%
%
%
%
%   Finally we also output a matlab script with the name OUT_NAME.m.
%
%   This is an M file generated in the working directory which automates 
%   the calling of save_graph, the binary and the load_graph functions.
%
%   The interface is:
%           [NEWVDATA, NEWADJ, NEWEDATA] = OUT_NAME(VDATA, ADJ, EDATA, ...
%                                                   OPTIONS, STRICT)
%
%      [VDATA, ADJ, EDATA]: graph representation. See 
%                           (OUT_NAME)_save_graph above for details
%
%      OPTIONS: A struct denoting the GraphLab parameters and the initial
%               schedule.
% 
%      OPTIONS.initial_schedule: (OUT_NAME)_save_graph above for details
%
%      OPTIONS.ncpus: A positive numeric denoting #threads to start. 
%
%      OPTIONS.scheduler: A string describing the scheduler and the 
%                         scheduler parameters. Identical to the 
%                         --scheduler=... parameter when running a 
%                         GraphLab program. 
%                         Run (OUT_NAME)_binary --help for details
%
%      OPTIONS.scope: A string describing the scope type. Should be 
%                     either 'edge', 'vertex' or 'full'. Identical to the 
%                         --scope=... parameter when running a 
%                         GraphLab program. 
%
%   The script first calls (OUT_NAME)_save_graph with the arguments VDATA,
%   ADJ, EDATA, OPTIONS.initial_schedule and STRUCT. The graph file name 
%   GRAPHFILE is automatically generated. The current directory must be 
%   writeable.  Next, it will invoke the program (OUT_NAME)_binary
%   passing it the rest of the parameters in OPTIONS. Finally the output
%   of the binary is read back using (OUT_NAME)_load_graph and returned.
%   
%
% Prereqs:
%   
%   This script assumes a *nix environment. It requires
%     -- MATLAB 2010b or later
%     -- g++ 4.3 or earlier if MATLAB 2010b.
%     -- make
%     -- A build/installation of GraphLab
%     -- If installed, this function must be located within the 
%        include/graphlab/matlab directory.
%
% For more details, see the README file. 

function compile_update_function(updatefunctions, exvertex_, exedge_, ...
                                 gllibdir, workingdirectory, outputmname, ...
                                 optimizationlevel, parcompiles, notcmalloc)

if (~exist('outputmname', 'var')) 
    outputmname = 'test';
end

if (~exist('optimizationlevel', 'var')) 
    optimizationlevel = 0;
end

if (~exist('parcompiles', 'var')) 
    parcompiles = 4;
end


if (~exist('notcmalloc', 'var')) 
    notcmalloc = 0;
end

if (~isscalar(parcompiles) || parcompiles < 1) 
    error('Parallel compiles must be >= 1');
    return;
end


if (~isscalar(notcmalloc)) 
    error('notcmalloc must be scalar');
    return;
end

[st,~] = dbstack('-completenames');
mfiledirectory = st.file;
% get the graphlab/matlab directory
slashes = strfind(mfiledirectory, '/');
glmatlabdir = mfiledirectory(1:(slashes(end)-1));
linkfndirectory = [glmatlabdir, '/link_functions'];
% get the graphlab directory
slashes = strfind(glmatlabdir, '/');
gldir = glmatlabdir(1:(slashes(end)-1));

% make the working directory a complete path
workingdirectory = [pwd '/' workingdirectory];

disp(['Graphlab source directory in :', gldir]);
disp(['Graphlab matlab source directory in :', glmatlabdir]);
disp(['Working Directory:', workingdirectory]);
% remember the current directory
olddir = pwd;
oldpath = path;

% ---------- Stage 0 ------------------------------------
% Workaround: There is an interesting issue if the edge
% structure is exactly equivalent to the vertex structure, or a substructure
% of the vertex structure.  (or vice versa)
% for instance:
%
% temp.a = []
% vdata.b = temp;
%
% edata = temp;
%
% essentially, when instantiating the 'b' term in the vertex structure, 
% emlc will complain that there are 2 candidate structures which match 
% the required definition. Either 'emx_edgedata' or 'some random string'
% 
% The solution is to purturb both structures with a random useless field
% which guarantees that they do not overlap definitions
%--------------------------------------------------------
fprintf('\n\n');
fprintf('-----------------------\n');
disp(['Stage 0: EMLC workarounds']);
fprintf('-----------------------\n');
if (isstruct(exvertex_))
    rstring = char(randi(26,[1,20])+'a'-1);
    disp(['adding field ' rstring ' to vertex struct']);
    exvertex_.(rstring) = 0.0;
end
if (isstruct(exedge_))
    rstring = char(randi(26,[1,20])+'a'-1);
    disp(['adding field ' rstring ' to edge struct']);
    exedge_.(rstring) = 0.0;
end


% ---------- Stage 1 ------------------------------------
% Translate update functions and interface functions to C
% -------------------------------------------------------
fprintf('\n\n');
fprintf('-------------------------------------------------------------\n');
disp(['Stage 1: Verify Datatype Limitations and Matlab Code Generation']);
fprintf('-------------------------------------------------------------\n');
[exvertex,status, genv] = gl_emx_typecheck(exvertex_, 'vdata');
if (status == 0)
    disp('Vertex data uses capabilities which exceed the GraphLab/EMLC specification');
    disp('Please see the documentation for details');
    return;
end

[exedge, status gene] = gl_emx_typecheck(exedge_, 'edata');
if (status == 0)
    disp('Vertex data uses capabilities which exceed the GraphLab/EMLC specification');
    disp('Please see the documentation for details');
    return;
end

disp('Checks passed.');
try

% ---------- Stage 2 --------------------
% EMLC Compilation
% ---------------------------------------
fprintf('\n\n');
fprintf('---------------------------\n');
disp(['Stage 2: EMLC Generation']);
fprintf('---------------------------\n');
% create and shift the the working directory
if (~isdir(workingdirectory))
    mkdir(workingdirectory);
end
if (~isdir(workingdirectory))
    disp('Failed to create working directory!');
    return;
end
cd(workingdirectory);

% make sure the old pwd is still in the path so that we can pick up the
% update functions
addpath(olddir);
% drop the link functions directory since we are making new ones
if (~isempty(strfind(path, linkfndirectory)))
    rmpath(linkfndirectory);
end
% keep the directory of this function in the path
addpath(glmatlabdir);

disp(['Generating Matlab<->Graphlab Link Functions']);

% generate the get vertexdata and get edge data functions
generate_link_functions(exvertex_, exedge_, genv, gene);

disp(['EMLC generation']);

% construct the example input for update functions
%function update_function(currentvertex,  % scalar integer
%                         inedges,        % array of in edge ids
%                         inv,            % array of sources vertices for each in edge id
%                         outedges,       % array of out edge ids
%                         outv,           % array of target vertices for each in edge id
%                         handle)         % scalar integer                    

exinput = ['-eg {uint32(0), ' ...               % current vertex
          'emlcoder.egs(uint32(0),[Inf]), ' ... % in edges
          'emlcoder.egs(uint32(0),[Inf]), ' ... % inV
          'emlcoder.egs(uint32(0),[Inf]), ' ... % out edges
          'emlcoder.egs(uint32(0),[Inf]), ' ... % outV
          'double(0.0)} '];                       % handle
          
cfg = emlcoder.CompilerOptions;
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.EnableVariableSizing = true;

emlcstring = ['emlc -s cfg -c -T RTW -d . -o updates -report'];
% append all the update functions
for i = 1:length(updatefunctions)
    emlcstring = [emlcstring ' ' updatefunctions{i} ' ' exinput];
end

% append the additional functions
emlcstring = [emlcstring ' datatype_identifier -eg {exvertex, exedge, emlcoder.egs(double(0), [Inf])}'];
emlcstring = [emlcstring ' get_vertex_data -eg {double(0), uint32(0)}'];
emlcstring = [emlcstring ' get_edge_data -eg {double(0), uint32(0)}'];
emlcstring = [emlcstring ' set_vertex_data -eg {double(0), uint32(0), exvertex}'];
emlcstring = [emlcstring ' set_edge_data -eg {double(0), uint32(0), exedge}'];
emlcstring = [emlcstring ' add_task -eg {double(0), uint32(0), emlcoder.egs(''a'', [Inf]), double(0)}'];
emlcstring = [emlcstring ' matlab_link.h'];
disp(['Issuing command: ' emlcstring]);
eval(emlcstring);
disp('EMLC Done');

% ---------- Stage 3 ---------------------
% Generate mxArray <-> emxArray converters
% ----------------------------------------
fprintf('\n\n');
fprintf('--------------------------------------------------------------\n');
disp('Stage 3: Generating Matlab <-> Embedded Matlab datatype converters');
fprintf('--------------------------------------------------------------\n');
cmd = ['python ' glmatlabdir '/mxarray_to_emlc.py updates_types.h ' glmatlabdir '/graphlab_options_struct.h > mx_emx_converters.hpp'];
disp(['Issuing command: ' cmd]);
[~,res] = system(cmd);
if (~isempty(res))
    disp 'Compilation Failed. ';
    error(res);
end

disp('Generating Update Function list:');
generate_update_function_list(updatefunctions);

% ---------- Stage 4 --------------------
% Compilation
% ---------------------------------------
fprintf('\n\n');
fprintf('---------------------\n');
disp('Stage 4: MEX compilation');
fprintf('---------------------\n');
% pick up all the generated c files
allcfiles = dir('*.c');
str = '';
for i = 1:length(allcfiles)
    str = [str allcfiles(i).name ' '];
    allcfilescellarray{i} = allcfiles(i).name;
end

largearraydims = [];
if (strcmp(computer, 'GLNXA64') ~= 1)
    largearraydims = ['-DMX_COMPAT_32'];
end

% 
% compilestring = ['mex -g '  ...
%   '-I. -I/usr/include ' ...
%  '-I' glmatlabdir ' '...
%  '-I' gldir '/../ ' ...
%  '-L' gllibdir ' ' ...
% '-L' gllibdir '/extern/metis/libmetis ' ...
% '-L' gllibdir '/extern/metis/GKlib ' ...
%  '-lgraphlab_pic ' ...
%  '-lgraphlab_util_pic ' ...
%  '-lgraphlab_metis_pic ' ...
%  '-lgraphlab_GKlib_pic ' ...
%  '-cxx -v ' ...
%  largearraydims ...
%  'CC="g++" ' ...
%  '-DBOOST_UBLAS_ENABLE_PROXY_SHORTCUTS ' ...
%  '-D_SCL_SECURE_NO_WARNINGS ' ...
%  '-D_CRT_SECURE_NO_WARNINGS ' ...
%  '-D_SECURE_SCL=0 ' ...
%  '-DMEX_COMPILE ' ...   
%  'CXXFLAGS="$CXXFLAGS -g -fpic -Wall -fno-strict-aliasing" ' ...
%  'CFLAGS="-D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread" ' ...
%  '-output test_mex ' ...
%  glmatlabdir '/graphlab_mex.cpp ' ...
%  glmatlabdir '/matlab_link.cpp ' ...
%  glmatlabdir '/graphlab_mex_parse.cpp ' ...
%  glmatlabdir '/graphlab_mex_output.cpp ' ...
%  glmatlabdir '/update_function_generator.cpp ' ...
%  str ];

%build up a list of all the files to compile
%compilefiles = allcfilescellarray;
compilefiles = {};
%compilefiles{length(compilefiles)+1} = [glmatlabdir '/matlab_link.cpp '];
compilefiles{length(compilefiles)+1} = [glmatlabdir '/mex_save_graph.cpp '];
compilefiles{length(compilefiles)+1} = [glmatlabdir '/graphlab_mex_parse.cpp '];
compilefiles{length(compilefiles)+1} = [glmatlabdir '/graphlab_mex_output.cpp '];
%compilefiles{length(compilefiles)+1} = [glmatlabdir '/update_function_generator.cpp '];

includepaths = {};
includepaths{length(includepaths)+1} = '.';
includepaths{length(includepaths)+1} = '/usr/include';
includepaths{length(includepaths)+1} = [gldir '/../'];
includepaths{length(includepaths)+1} = glmatlabdir;

libpaths = {};
libpaths{length(libpaths)+1} = gllibdir ;
if (exist([gllibdir '/extern/metis/libmetis'],'dir'))
    libpaths{length(libpaths)+1} = [gllibdir '/extern/metis/libmetis'] ;
end
if (exist([gllibdir '/extern/metis/GKlib'],'dir'))
    libpaths{length(libpaths)+1} = [gllibdir '/extern/metis/GKlib'] ;
end

optlevelstring = '';
if optimizationlevel == 1
    optlevelstring = '-O1';
elseif optimizationlevel == 2
    optlevelstring = '-O2';
elseif optimizationlevel >= 3
    optlevelstring = '-O3';
end
generate_makefile( [workingdirectory '/Makefile_save'], ...
                       [outputmname '_save_graph'], ...
                       compilefiles, ...
                       includepaths, ...
                       libpaths, ...
                       ['-DMEX_COMPILE -g -Wall -fno-strict-aliasing -D_GNU_SOURCE -fexceptions -fno-omit-frame-pointer -pthread ' optlevelstring ' ' largearraydims], ...
                       '-lgraphlab_pic -lgraphlab_util_pic -lgraphlab_metis_pic -lgraphlab_GKlib_pic', ...
                       true);

compilefiles = {};
%compilefiles{length(compilefiles)+1} = [glmatlabdir '/matlab_link.cpp '];
compilefiles{length(compilefiles)+1} = [glmatlabdir '/mex_load_graph.cpp '];
compilefiles{length(compilefiles)+1} = [glmatlabdir '/graphlab_mex_parse.cpp '];
compilefiles{length(compilefiles)+1} = [glmatlabdir '/graphlab_mex_output.cpp '];
%compilefiles{length(compilefiles)+1} = [glmatlabdir '/update_function_generator.cpp '];

generate_makefile( [workingdirectory '/Makefile_load'], ...
                       [outputmname '_load_graph'], ...
                       compilefiles, ...
                       includepaths, ...
                       libpaths, ...
                       ['-DMEX_COMPILE -g -Wall -fno-strict-aliasing -D_GNU_SOURCE -fexceptions -pthread ' optlevelstring ' ' largearraydims], ...
                       '-lgraphlab_pic -lgraphlab_util_pic -lgraphlab_metis_pic -lgraphlab_GKlib_pic', ...
                       true);

% check for tcmalloc
tcmalloclib = '';
if (notcmalloc == 0)
    [~,r] = system('whereis -b libtcmalloc.a');
    if (~isempty(strfind(r,'libtcmalloc.a')))
        tcmalloclib = '-ltcmalloc';
    else
        disp('tcmalloc not found (using whereis). It should be installed to obtain optimal performance.');
    end
    
end
   

bincompilefiles = allcfilescellarray;
bincompilefiles{length(bincompilefiles)+1} = [glmatlabdir '/binary_stage.cpp '];
bincompilefiles{length(bincompilefiles)+1} = [glmatlabdir '/matlab_link.cpp '];
bincompilefiles{length(bincompilefiles)+1} = [glmatlabdir '/update_function_generator.cpp '];

generate_makefile( [workingdirectory '/Makefile_bin'], ...
                       [outputmname '_binary'], ...
                       bincompilefiles, ...
                       includepaths, ...
                       libpaths, ...
                       ['-g -Wall -fno-strict-aliasing -fno-omit-frame-pointer -pthread ' optlevelstring ' ' largearraydims], ...
                       ['-lgraphlab -lgraphlab_util -lgraphlab_metis -lgraphlab_GKlib -lpthread -lboost_program_options ' tcmalloclib], ...
                       false, ...
                       'b_');

%unix(['make -j' num2str(parcompiles)]);
                       
%disp(['Issuing command: ' compilestring]);
%eval(compilestring);

catch exception
    % revert 
    path(oldpath);
    cd(olddir);
    disp(getReport(exception,'extended'));
    rethrow(exception);
end
fprintf('\n\n');
disp(['Complete.']);

% revert 
path(oldpath);
cd(olddir);
rehash;
end