function test_update_function(currentvertex, inedges, inv, outedges, outv, handle)   
    vdata = get_vertex_data(handle, currentvertex);
    vdata.b = 'moo';
    temp = zeros(4);
    eml.varsize('temp');
    for i = 1:length(inedges)
        edata = get_edge_data(handle, inedges(i));
        temp = temp + double(edata.b);
        edata.b = uint16(zeros(4));
        set_edge_data(handle, inedges(i), edata);
    end
    vdata.c(1).a = temp;
    set_vertex_data(handle, currentvertex, vdata);
end
