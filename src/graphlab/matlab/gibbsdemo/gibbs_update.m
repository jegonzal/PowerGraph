function gibbs_update(currentvertex, inedges, inv, outedges, outv, handle) %#eml
    % get the current vertex data
    vdata = get_vertex_data(handle, currentvertex);
    
    % base probability
    logp = vdata.logunary;
    for i = 1:length(inedges)
        adjacent_vertex = get_vertex_data(handle, inv(i));
        edge_potential = get_edge_data(handle,inedges(i));
        logp = logp + edge_potential(adjacent_vertex.sample,:);
    end
    % do some initial normalization so we dont get over/underflows 
    % after taking exp
    logp = logp - min(logp);
    p = exp(logp);
    % draw the sample
    vdata.sample = rand_multinomial(p);
    % set the new vertex data
    set_vertex_data(handle, currentvertex, vdata); 
end
