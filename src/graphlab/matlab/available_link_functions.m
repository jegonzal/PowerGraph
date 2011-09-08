function linkfn = available_link_functions
linkfn = [];
% A list of all the available link functions.
% each of this functions must have a matching [name].template file
% which issues the call to a C function

linkfn(1).name = 'datatype_identifier';
linkfn(1).args = 'exvertex, exedge, emlcoder.egs(handle_type__)';

linkfn(2).name ='get_vertex_data';
linkfn(2).args = 'emlcoder.egs(handle_type__), uint32(0)';

linkfn(3).name = 'get_edge_data';
linkfn(3).args = 'emlcoder.egs(handle_type__), uint32(0)';

linkfn(4).name = 'set_vertex_data';
linkfn(4).args = 'emlcoder.egs(handle_type__), uint32(0), exvertex';

linkfn(5).name = 'set_edge_data';
linkfn(5).args = 'emlcoder.egs(handle_type__), uint32(0), exedge';

linkfn(6).name = 'add_task';
linkfn(6).args = 'emlcoder.egs(handle_type__), uint32(0), emlcoder.egs(''a'', [Inf]), double(0)';


linkfn(7).name = 'rand_bernoulli';
linkfn(7).args = 'double(0)';

%linkfn(8).name = 'rand_bernoulli_fast';
%linkfn(8).args = 'double(0)';

%linkfn(9).name = 'rand_double';
%linkfn(9).args = '';

%linkfn(10).name = 'rand_gamma';
%linkfn(10).args = 'double(0)';
% 
% linkfn(11).name = 'rand_gaussian';
% linkfn(11).args = 'double(0), double(0)';
% 
% linkfn(12).name = 'rand_int';
% linkfn(12).args = '';
% 
% linkfn(13).name = 'rand_int_uniform';
% linkfn(13).args = 'uint32(0)';
% 
% linkfn(14).name = 'rand_int_uniform_fast';
% linkfn(14).args = 'uint32(0)';
% 
% linkfn(15).name = 'rand_multinomial';
% linkfn(15).args = 'emlcoder.egs(double(0.0), [Inf])';