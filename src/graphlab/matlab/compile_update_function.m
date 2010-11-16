function compile_update_function(updatefunctions, exvertex_, exedge_, gllibdir, workingdirectory)
[st,~] = dbstack('-completenames');
mfiledirectory = st.file;
% get the graphlab/matlab directory
slashes = strfind(mfiledirectory, '/');
glmatlabdir = mfiledirectory(1:(slashes(end)-1));
linkfndirectory = [glmatlabdir, '/link_functions'];
% get the graphlab directory
slashes = strfind(glmatlabdir, '/');
gldir = glmatlabdir(1:(slashes(end)-1));


disp(['Graphlab source directory in :', gldir]);
disp(['Graphlab matlab source directory in :', glmatlabdir]);
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
rmpath(linkfndirectory);
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
emlcstring = [emlcstring ' datatype_identifier -eg {exvertex, exedge}'];
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
cmd = ['python ' glmatlabdir '/mxarray_to_emlc.py updates_types.h > mx_emx_converters.hpp'];
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
end

largearraydims = [];
if (computer == 'GLNXA64')
    largearraydims = ['-largeArrayDims '];
end
compilestring = ['mex -g '  ...
  '-I. -I/usr/include ' ...
 '-I' glmatlabdir ' '...
 '-I' gldir '/../ ' ...
 '-L' gllibdir ' ' ...
'-L' gllibdir '/extern/metis/libmetis ' ...
'-L' gllibdir '/extern/metis/GKlib ' ...
 '-lgraphlab_pic ' ...
 '-lgraphlab_util_pic ' ...
 '-lgraphlab_metis_pic ' ...
 '-lgraphlab_GKlib_pic ' ...
 '-cxx -v ' ...
 largearraydims ...
 'CC="g++" ' ...
 '-DBOOST_UBLAS_ENABLE_PROXY_SHORTCUTS ' ...
 '-D_SCL_SECURE_NO_WARNINGS ' ...
 '-D_CRT_SECURE_NO_WARNINGS ' ...
 '-D_SECURE_SCL=0 ' ...
 '-DMEX_COMPILE ' ...   
 'CXXFLAGS="$CXXFLAGS -g -fpic -Wall -fno-strict-aliasing" ' ...
 'CFLAGS="-D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread" ' ...
 '-output test_mex ' ...
 glmatlabdir '/graphlab_mex.cpp ' ...
 glmatlabdir '/matlab_link.cpp ' ...
 glmatlabdir '/graphlab_mex_parse.cpp ' ...
 glmatlabdir '/graphlab_mex_output.cpp ' ...
 glmatlabdir '/update_function_generator.cpp ' ...
 str ];


disp(['Issuing command: ' compilestring]);
eval(compilestring);

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
end