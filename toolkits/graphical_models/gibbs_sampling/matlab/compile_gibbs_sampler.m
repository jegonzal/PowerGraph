%% Compile all mex files
% This script compiles the mex files needed to run the parallel sampling
% algorithms.

graphlab_path='../../..';
graphlab_bin_path=[graphlab_path, '/release'];
pgibbs_bin_path=[graphlab_bin_path, '/demoapps/pgibbs'];
graphlab_include_path = [graphlab_path, '/src'];
graphlab_link_path = [graphlab_bin_path, '/src/graphlab'];
cxx_flags = ['CXXFLAGS=', ...
   '"-fPIC -Wall -O3 -pthread -fexceptions -fno-omit-frame-pointer ', ...
   '-fopenmp"'];

 
%% If release folder does not exist run configure
if(~exist(graphlab_bin_path, 'dir')) 
  disp('Configure was not yet run running config now.');
  [errorstatus, result] = ...
    system(['cd ', graphlab_path, ';', ' ./configure'])
  if(errorstatus) 
    error('Error running config!');
  end
end

%% Compile the pgibbs library needed for the mex file
[errorstatus, result] = ...
  system(['cd ', pgibbs_bin_path, ';', ' make -j2'])
if(errorstatus) 
  error('Error compiling pgibbs!');
end

  
%% Do the compilation
compiler_type_flags = '';
if(ismac()) 
  disp('We require gcc 4.2 on mac');
  compiler_type_flags = 'LD=gcc-4.2 CC=gcc-4.2 CXX=g++-4.2';
end 

compilestr = ...
   ['mex ', compiler_type_flags, ' ', ...
    '-largeArrayDims', ' ', ...
    cxx_flags, ' ', ...
    'gibbs_sampler_impl.cpp', ' ', ...
    '-I', graphlab_include_path, ' ', ...
    '-L', graphlab_link_path, ' ', ...
    '-L', graphlab_link_path, '/extern/metis/GKlib', ' ', ...
    '-L', graphlab_link_path, '/extern/metis/libmetis', ' ', ...
    '-L', pgibbs_bin_path, ' ', ...
    '-lpgibbs_pic', ' ', ...
    '-lgraphlab_pic ', ' ', ...
    '-lgomp'];

disp(compilestr);
eval(compilestr);
disp('Finished!');
