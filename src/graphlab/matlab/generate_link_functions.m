% generate the functions the update functions can call.
% get_vertex_data, set_vertex_data, get_edge_data, set_edge_data
% as well as the datatype_identifier
% exvertex_: example vertex data
% exedge_: example edge data
% genv: matlab code that generates exvertex_ as well as specify varsize for
%       all arrays in exvertex_
% gene: matlab code that generates exedge_ as well as specify varsize for
%       all arrays in exedge_
function generate_link_functions(exvertex_, exedge_, genv, gene)

% get_vertex_data can be called by the user function to get the data
% on a vertex. This call works by redirecting to a C function
% emx_get_vertex_data
f = fopen(['get_vertex_data.m'], 'w');
fprintf(f, 'function vdata = get_vertex_data(handle, vertex) %%#eml\n');
fprintf(f, [genv, '\n']);
if (isstruct(exvertex_))
    fprintf(f, 'eml.cstructname(vdata, ''emx_vertexdata'', ''extern'');\n');
end
fprintf(f, 'eml.ceval(''emx_get_vertex_data'', handle, vertex, eml.ref(vdata));\n');
fprintf(f, 'end\n');
fclose(f);

% get_edge_data can be called by the user function to get the data
% on an edge. This call works by redirecting to a C function
% emx_get_edge_data
f = fopen(['get_edge_data.m'], 'w');
fprintf(f, 'function edata = get_edge_data(handle, edge) %%#eml\n');
fprintf(f, [gene, '\n']);
if (isstruct(exedge_))
    fprintf(f, 'eml.cstructname(edata, ''emx_edgedata'', ''extern'');\n');
end
fprintf(f, 'eml.ceval(''emx_get_edge_data'', handle, edge, eml.ref(edata));\n');
fprintf(f, 'end\n');
fclose(f);

% datatype_identifier is an empty function that allows the C++
% side to identify the vertex and edge data types
f = fopen(['datatype_identifier.m'], 'w');
fprintf(f, 'function datatype_identifier(vdata, edata) %%#eml\n');
if (isstruct(exvertex_))
    fprintf(f, 'eml.cstructname(vdata, ''emx_vertexdata'');\n');
end
if (isstruct(exedge_))
    fprintf(f, 'eml.cstructname(edata, ''emx_edgedata'');\n');
end
fprintf(f, 'end\n');
fclose(f);



% set_vertex_data can be called by the user function to set the data
% on a vertex. This call works by redirecting to a C function
% emx_set_vertex_data
f = fopen(['set_vertex_data.m'], 'w');
fprintf(f, 'function set_vertex_data(handle, vertex, vdata) %%#eml\n');
if (isstruct(exvertex_))
    fprintf(f, 'eml.cstructname(vdata, ''emx_vertexdata'', ''extern'');\n');
end
fprintf(f, 'eml.ceval(''emx_set_vertex_data'', handle, vertex, eml.ref(vdata));\n');
fprintf(f, 'end\n');
fclose(f);

% set_edge_data can be called by the user function to get the data
% on an edge. This call works by redirecting to a C function
% emx_set_edge_data
f = fopen(['set_edge_data.m'], 'w');
fprintf(f, 'function set_edge_data(handle, edge, edata) %%#eml\n');
if (isstruct(exedge_))
    fprintf(f, 'eml.cstructname(edata, ''emx_edgedata'', ''extern'');\n');
end
fprintf(f, 'eml.ceval(''emx_set_edge_data'', handle, edge, eml.ref(edata));\n');
fprintf(f, 'end\n');
fclose(f);


% add task
f = fopen(['add_task.m'], 'w');
fprintf(f, 'function add_task(handle, vertex, functionname, priority) %%#eml\n');
fprintf(f, '  eml.ceval(''emx_add_task'', handle, vertex, [functionname, 0], priority)\n');
fprintf(f, 'end');
fclose(f);


end
