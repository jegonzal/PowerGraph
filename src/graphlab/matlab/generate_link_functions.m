% generate the functions the update functions can call.
% get_vertex_data, set_vertex_data, get_edge_data, set_edge_data
% as well as the datatype_identifier
% exvertex_: example vertex data
% exedge_: example edge data
% genv: matlab code that generates exvertex_ as well as specify varsize for
%       all arrays in exvertex_
% gene: matlab code that generates exedge_ as well as specify varsize for
%       all arrays in exedge_
function generate_link_functions(exvertex_, exedge_, genv, gene, templatedirectory)
substs.GENVERTEX = genv;
substs.GENEDGE = gene;
if (isstruct(exvertex_))
     substs.EXTERN_VERTEX_STRUCT = 'eml.cstructname(vdata, ''emx_vertexdata'', ''extern'');';
     substs.VERTEX_STRUCT = 'eml.cstructname(vdata, ''emx_vertexdata'');';
else
     substs.EXTERN_VERTEX_STRUCT = '';
     substs.VERTEX_STRUCT ='';
end
 
if (isstruct(exedge_))
     substs.EXTERN_EDGE_STRUCT = 'eml.cstructname(edata, ''emx_edgedata'', ''extern'');';
     substs.EDGE_STRUCT = 'eml.cstructname(edata, ''emx_edgedata'');';
else
     substs.EXTERN_EDGE_STRUCT = '';
     substs.EDGE_STRUCT = '';
end

substs.HANDLE_STRUCT = 'eml.cstructname(handle, ''HANDLE_TYPE'');';

file_template_substitution([templatedirectory, '/get_vertex_data.template'], ...
                           'get_vertex_data.m', ...
                           substs);

file_template_substitution([templatedirectory, '/get_edge_data.template'], ...
                           'get_edge_data.m', ...
                           substs);

file_template_substitution([templatedirectory, '/set_vertex_data.template'], ...
                           'set_vertex_data.m', ...
                           substs);

file_template_substitution([templatedirectory, '/set_edge_data.template'], ...
                           'set_edge_data.m', ...
                           substs);

file_template_substitution([templatedirectory, '/datatype_identifier.template'], ...
                           'datatype_identifier.m', ...
                           substs);

file_template_substitution([templatedirectory, '/add_task.template'], ...
                           'add_task.m', ...
                           substs);
end
