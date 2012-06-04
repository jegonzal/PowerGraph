%% Parallel Gibbs sampler
% The parallel gibbs sampler is an optimized a c++ implementation of
% the discrete Gibbs samplers which uses multiple threads to
% accelerate the generation of a single sampling chain.  The parallel
% Gibbs sampler implements two algorithms described in the paper:
%
%   Parallel Gibbs Sampling: From Colored Fields to Think Junction Trees
%     by Joseph Gonzalez, Yucheng Low, Arthur Gretton, and Carlos Guestrin
%
% The first algorithm is the Chromatic sampler which is a direct
% parallelization of the classic Gibbs sampler.  The second algorithm
% is the Splash Gibbs sampler which incrementally builds thin junction
% trees.
%
% To use this function you must first construct a discrete factor
% graph which is simply a cell array of table factors:
%
%   factor{1} = table_factor( [1,2], log(rand(3,4)) );
%   factor{2} = table_factor( [2,3], log(rand(4,2)) );
%
% This creates a factorized model (with random tables) over the
% variables 1, 2, and 3. We can then run the CHROMATIC sampler by
% calling:
%
%   options.alg_type = 'CHROMATIC';
%   options.nsamples = 100;
%   options.nskip    = 10;
%   [samples, nupdates, nchanges, marginals] = ...
%      gibbs_sampler(factors, options); 
%
%
% 
% Arguments:
%   factors: a cell array of factors constructed using the table_factor 
%     function.
%   options: a struct with the following fields:
%    * alg_type: [Default: 'CHROMATIC'] A string either 'CHROMATIC' or
%        'SPLASH'.  For relatively fast mixing models the 'CHORMATIC'
%        algorithm is simpler and faster.  For slowly mixing models
%        use the 'SPLASH' algorithm.  In this case additional options
%        will need to be set.
%    * nsamples: [Default: 10] The number of joint samples to collect.
%    * nskip: [Default: 10] The number of samples to skip between
%        joint samples.  Because of the asynchronous nature of the
%        algorithms more or than nskip samples may actually be skipped
%        in practice.  In the 'SPLASH' algorithm nskip * nvariables
%        single variable updates are computed before the next joint
%        sample is constructed.
%    * ncpus: [Default: 2] The number of threads to use when running
%        the inference algorithm. The number of cpus should be less
%        than the number of variables and ideally not much larger than
%        the number of processors.
%    * treewidth: [Default: 3] The treewidth of the junction trees
%        constructed using the Splash sampler.
%    * treeheight: [Default: maxint] The largest height of a tree
%    * treesize: [Default: maxint] The largest possible size of a
%        tree
%    * priorities: [Default: false] Use priorities when
%       constructing the splash trees
%    * checkargs: [Default: True] Determines if the arguments are
%        checked before calling the C++ code.  While we do additional
%        argument checking withing the C++ code it is often easier to
%        debug broken factors from within the matlab code.  However
%        for the fastest performance disable checkargs (set to false).
%
% Return Arguments:
%   samples: nvars * nsamples matrix of joint assingments
%   nupdates: nvars * nsamples number of times each variable was updated.
%   nchanges: nvars * nsamples the number of times the variable's 
%      assignment changed values
%   beliefs: nvars * nsamples cell array of vectors represent the 
%      Rao-Blackwellized marginal estimates for each variable.
%
% See Also: table_factor
%
% This actual c++ mex function is provided in gibbs_sampler_impl.cpp
% which can be compiled by running compile_gibbs_sampler.m.
%
function [samples, nupdates, nchanges, marginals] = ...
      gibbs_sampler(factors, options)    

  %% Check the arguments
  if(~iscell(factors))
    error('The factors argument must be a cell array of table_factors');
  end
  % Define default options
  if(~exist('options', 'var'))
    options.alg_type = 'CHROMATIC';
  end
  if(~isfield(options, 'alg_type'))
    options.alg_type = 'CHROMATIC';
  end
  if(~isfield(options, 'nsamples'))
    options.nsamples = 10;
  end
  if(~isfield(options, 'nskip'))
    options.nskip = 10;
  end
  if(~isfield(options, 'ncpus'))
    options.ncpus = 2;
  end
  if(~isfield(options, 'treewidth'))
    options.treewidth = 3;
  end
  if(~isfield(options, 'treeheight'))
    options.treeheight = double(intmax());
  end
  if(~isfield(options, 'treesize'))
    options.treesize = double(intmax());
  end
  if(~isfield(options, 'priorities'))
    options.priorities = false;
  end
  if(~isfield(options, 'vanish'))
    options.vanish = 10;
  end
  if(~isfield(options, 'checkargs'))
    options.checkargs = true;
  end
  if(~isfield(options, 'save_alchemy'))
    options.save_alchemy = false;
  end

  options.treewidth = double(options.treewidth);
  if(options.treewidth > 32) 
     error('Treewidth must be less than 32');
  end
  if(options.treewidth < 1) 
     error('Treewidth must be at least 1');
  end

  if(options.checkargs) 
    max_var = 0;
    %% Check the factors data structure
    for i = 1:length(factors)
      if(~isfield(factors{i}, 'vars'))
        disp(factors{i});
        error('Factor %d does not contain the field vars', i);
      end
      if(~strcmp(class(factors{i}.vars), 'uint32'))
        disp(factors{i});
        error(['Factor ', num2str(i), ...
               ' has variables of type ', ...
               class(factors{i}.vars), ...
               ' when they should be of type uint32.']);    
      end
      if(~isfield(factors{i}, 'logP'))
        disp(factors{i});
        error('Factor %d does not contain the field logP', i);
      end
      if(~strcmp(class(factors{i}.logP), 'double'))
        disp(factors{i});
        error(['Factor ', num2str(i), ...
               ' has logP of type ', ...
               class(factors{i}.logP), ...
               ' when they should be of type double.']);    
      end    
      % Get the maximum variables
      max_var = max(max(factors{i}.vars(:)), max_var);
      if(min(factors{i}.vars) <= 0) 
        disp(factors{i});
        error('Factor %d has 0 valued variables', i);
      end   
    end

    %% check all the variables have consistent sizes;
    vars = 1:max_var;
    var_sizes = zeros(max_var, 1);
    for i = 1:length(factors)
      current_sizes = var_sizes(factors{i}.vars);
      dims = size(factors{i}.logP)';
      dims = dims(dims > 1);
      if(length(dims) ~= length(current_sizes))
        error(['The number of variables %d in factor %d does not match ' ...
               'the number of dimensions %d.'], ...
              length(current_sizes), i, length(dims));
      end   
      ind = current_sizes > 0;
      % all elements in ind have been set and should just match
      if( ~isempty(find(current_sizes(ind) ~= dims(ind), 1)) ) 
        errorind = find(current_sizes(ind) ~= dims(ind), 1);
        error(['Variable %d has already been seen having size %d ', ...
               'but was just now observed to have size %d in factor %d.'], ...
              factors{i}.variables(ind(errorind)), ...
              current_sizes(ind(errorind)), ...
              dims(ind(errorind)), ...
              i);
      end
      var_sizes(factors{i}.vars(~ind)) = dims(~ind);
    end
    unset_vars = find(var_sizes(:) == 0);
    if(~isempty(unset_vars)) 
      error(['The following variables were not set correctly: ', ...
             mat2str(unset_vars)]);   
    end
  end
  
  %% Call the sampler
  if(nargout() <= 1)
    samples = gibbs_sampler_impl(factors, options);
  elseif(nargout() == 2)
    [samples, nupdates] = gibbs_sampler_impl(factors, options);
  elseif(nargout() == 3)
    [samples, nupdates, nchanges] = ...
      gibbs_sampler_impl(factors, options);
  elseif(nargout() == 4)
    [samples, nupdates, nchanges, marginals] = ...
      gibbs_sampler_impl(factors, options);
  end

end
    
    
    
    