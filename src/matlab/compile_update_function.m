function compile_update_function(updatefunctions, exvertex_, exedge_, workingdirectory)
% remember the current directory
olddir = pwd;
oldpath = path;

% ---------- Stage 0 ------------------------------------
% Bug check: There is an interesting issue if the edge
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
exvertex_.(char(randi(26,[1,20])+'a'-1)) = 0.0;
exedge_.(char(randi(26,[1,20])+'a'-1)) = 0.0;


% ---------- Stage 1 ------------------------------------
% Translate update functions and interface functions to C
% -------------------------------------------------------

[exvertex,status, genv] = gl_emx_typecheck(exvertex_, 'vdata');
if (status == 0)
    disp('Vertex data uses capabilities which exceed the GraphLab/EMLC specification');
    disp('Please see the documentation for details');
    return
end

[exedge, status gene] = gl_emx_typecheck(exedge_, 'edata');
if (status == 0)
    disp('Vertex data uses capabilities which exceed the GraphLab/EMLC specification');
    disp('Please see the documentation for details');
    return
end

try

% create and shift the the working directory
mkdir(workingdirectory);
cd(workingdirectory);
% make sure the old pwd is still in the path so that we can pick up the
% update functions
path(path, olddir);



% generate the get vertexdata and get edge data functions
generate_link_functions(exvertex_, exedge_, genv, gene);



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
          'uint32(0)} '];                       % handle
          
cfg = emlcoder.CompilerOptions;
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.EnableVariableSizing = true;

emlcstring = ['emlc -s cfg -c -T RTW -d . -o updates'];
% append all the update functions
for i = 1:length(updatefunctions)
    emlcstring = [emlcstring ' ' updatefunctions{i} ' ' exinput];
end

% append the additional functions
emlcstring = [emlcstring ' datatype_identifier -eg {exvertex, exedge}'];
emlcstring = [emlcstring ' get_vertex_data -eg {uint32(0), uint32(0)}'];
emlcstring = [emlcstring ' get_edge_data -eg {uint32(0), uint32(0)}'];
emlcstring = [emlcstring ' set_vertex_data -eg {uint32(0), uint32(0), exvertex}'];
emlcstring = [emlcstring ' set_edge_data -eg {uint32(0), uint32(0), exedge}'];
emlcstring = [emlcstring ' matlab_link.h'];
disp(emlcstring);
eval(emlcstring);


% ---------- Stage 2 --------------------
% Generate mxArray <-> emxArray converters
% ---------------------------------------
[unused,res] = system(['python ' olddir '/mxarray_to_emlc.py updates_types.h > generator.hpp']);
if (~isempty(res))
    disp 'Compilation Failed. ';
    error(res);
end

% ---------- Stage 3 --------------------
% Compilation
% ---------------------------------------
% pick up all the generated c files
allcfiles = dir('*.c');
str = '';
for i = 1:length(allcfiles)
    str = [str allcfiles(i).name ' '];
end

compilestring = ['mex -g '  ...
  '-I. -I/usr/include ' ...
 '-I' olddir ' '...
 '-cxx ' ...
 'CC="g++" ' ...
 '-DBOOST_UBLAS_ENABLE_PROXY_SHORTCUTS ' ...
 '-D_SCL_SECURE_NO_WARNINGS ' ...
 '-D_CRT_SECURE_NO_WARNINGS ' ...
 '-D_SECURE_SCL=0 ' ...
 '-DMEX_COMPILE ' ...   
 'CXXFLAGS="$CXXFLAGS -g -fpic" ' ...
 '-output test_mex ' ...
 olddir '/test_mex.cpp ' ...
 olddir '/matlab_link.cpp ' ...
 str ];

disp(compilestring);
eval(compilestring);

catch exception
    % revert 
    path(oldpath);
    cd(olddir);
    disp(getReport(exception,'extended'));
    rethrow(exception);
end



% revert 
path(oldpath);
cd(olddir);
end