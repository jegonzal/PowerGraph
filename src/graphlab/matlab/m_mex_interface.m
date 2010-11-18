% this function provides the general schema for calling a MatLab GraphLab
% program.
% vertexdata: cell array of vertex data
% adj_mat: (sparse) adjacency matrix where adj_mat[i][j] is an edge from vertex
%          i to vertex j and the data on the edge is edgedata(adjmat[i][j])
% edgedata: cell array of edge data
% options: is a struct describing start up options as well as the initial
%          schedule. The details of the options are provided below.
% Returns new graph  data on exit
%
%
%
% options.scheduler is the scheduler string. This can be for instance
%                    'sweep', or 'fifo' or 'colored(10)'.
%
% options.ncpus is the number of threads to start
%
% options.scope is the scope model. This should be either 'vertex','edge'
%               or 'full'
%
% options.initial_schedule describes the initial scheduling.
% options.initial schedule is a struct array. where each struct is of the
% form.
% 
% sched.update_function : the name of the update function as a string
% sched.vertices        : an array of vertices to update using this update
%                         function
% sched.priorities     :  a double array of of the same size as
%                         sched.vertices defining the priority for each
%                         vertex. Each value should be strictly > 0
%
% options.initial_schedule is a struct array of 'sched', allowing
% multiple update functions to be defined in the schedule

function [vertexdata, adj_mat, edgedata] = m_mex_interface(vertexdata, adj_mat, edgedata, options)
% cast vertices as uint32
for i = 1:length(options.initial_schedule)
    options.initial_schedule(i).vertices = uint32(options.initial_schedule(i).vertices);
end
graphfile = 'graph.tmp';
ret = mex_save_graph(vertexdata, adj_mat, edgedata, options, graphfile, 1);
if ret == false
    return
end

% construct the binary call
command = ['./binary --scheduler="' options.scheduler '" --ncpus=' num2str(options.ncpus) ...
           ' --scope=' options.scope ' --ingraphfile=graph.tmp --outgraphfile=graphout.tmp'];

unix(command);

[vertexdata, adj_mat, edgedata] = mex_load_graph(graphoutfile);
end