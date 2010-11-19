function bp_update(currentvertex, inedges, inv, outedges, outv, handle) %#eml
    vdata = get_vertex_data(handle, currentvertex);
    % update belief
    % compute new belief
    vdata.belief = vdata.unary;
    for i = 1:length(inedges)
        inedata = get_edge_data(handle, inedges(i));
        vdata.belief = vdata.belief .* inedata.msg;
        vdata.belief = vdata.belief / sum(vdata.belief);
    end
    set_vertex_data(handle, currentvertex, vdata);    
    % write out messages
    for i = 1:length(inedges)
        inedata = get_edge_data(handle, inedges(i));
        outedata = get_edge_data(handle, outedges(i));
        % get the out going message
        oldoutedatamsg = outedata.msg;
        outedata.msg = vdata.belief ./ inedata.msg;
        outedata.msg = outedata.msg / sum(outedata.msg);
        
        % outedata.msg is a row vector
        % multiply by the edge factor and normalize
        outedata.msg = outedata.msg * outedata.binary;
        outedata.msg = outedata.msg / sum(outedata.msg);
        set_edge_data(handle,outedges(i),outedata);
        % compute the residual;
        residual = sum(abs(outedata.msg - oldoutedatamsg));
        if (residual > 1E-5)
            add_task(handle, outv(i), 'bp_update', residual);
        end
    end
end
