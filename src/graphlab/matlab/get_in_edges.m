function arr = get_in_edges(handle, vertex) %#eml
    arr = int32([]);
    eml.varsize('arr');
    eml.ceval('get_in_edges_impl', handle, vertex, eml.wref(arr));
end