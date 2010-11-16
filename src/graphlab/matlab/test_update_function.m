function test_update_function(currentvertex, inedges, inv, outedges, outv, handle)   
    vdata = get_vertex_data(handle, currentvertex);
    vdata.b = 'moo';
    set_vertex_data(handle, currentvertex, vdata);
end
