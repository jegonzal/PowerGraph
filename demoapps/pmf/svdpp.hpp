#ifndef __SVD_HPP
#define __SVD_HPP

void calc_user_moviebag(vertex_data & vdata, edge_list & outs){

#ifdef GL_SVD_PP
    vdata.weight = zeros(D);
    foreach(graphlab::edge_id_t oedgeid, outs) {
       //edge_data & edge = scope.edge_data(oedgeid);
       const vertex_data  & pdata = scope.const_neighbor_vertex_data(scope.target(oedgeid)); 
       vdata.weight += pdata.weight;
    }
#endif //GL_SVD_PP
} 


#endif //__SVD_HPP
