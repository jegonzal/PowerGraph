% generate the link functions the update functions can call.
% exvertex_: example vertex data
% exedge_: example edge data
% genv: matlab code that generates exvertex_ as well as specify varsize for
%       all arrays in exvertex_
% gene: matlab code that generates exedge_ as well as specify varsize for
%       all arrays in exedge_
function generate_link_functions(exvertex_, exedge_, genv, gene, templatedirectory, link_functions)
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

for linkfn = link_functions
file_template_substitution([templatedirectory, '/' linkfn.name '.template'], ...
                           [linkfn.name '.m'], ...
                           substs);
end

end
