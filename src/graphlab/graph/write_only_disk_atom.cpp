#include <graphlab/graph/write_only_disk_atom.hpp>
#include <graphlab/graph/disk_atom.hpp>
#include <graphlab/graph/memory_atom.hpp>

namespace graphlab {

void write_only_disk_atom::play_back(graph_atom* atom) {
  f.close();
  std::ifstream fin(filename.c_str(), std::ios::binary);
  // flush the commands
  iarchive iarc(fin);
  
  while(fin.good()) {
    char command;
    fin >> command;
    if (fin.fail()) break;
    if (command == 'a') {
      vertex_id_type vid; uint16_t owner;
      iarc >> vid >> owner;
      atom->add_vertex_skip(vid, owner);
    } else if (command == 'b') {
      vertex_id_type vid; uint16_t owner;
      iarc >> vid >> owner;
      atom->add_vertex_skip(vid, owner);
    } else if (command == 'c') {
      vertex_id_type vid; uint16_t owner; std::string data;
      iarc >> vid >> owner >> data;
      atom->add_vertex_with_data(vid, owner, data);
    } else if (command == 'f') {
      vertex_id_type src; vertex_id_type target;
      std::string data;
      iarc >> src >> target >> data;
      atom->add_edge_with_data(src, target, data);
    } else if (command == 'g') {
      vertex_id_type vid; uint16_t owner;
      iarc >> vid >> owner;
      atom->set_vertex(vid, owner);
    } else if (command == 'h') {
      vertex_id_type vid; uint16_t owner; std::string data;
      iarc >> vid >> owner >> data;
      atom->set_vertex_with_data(vid, owner, data);
    } else if (command == 'j') {
      vertex_id_type src; vertex_id_type target; std::string data;
      iarc >> src >> target >> data;
      atom->set_edge_with_data(src, target, data);
    } else if (command == 'd') {
      vertex_id_type src; vertex_id_type target; std::string data;
      uint16_t srcowner, targetowner;
      iarc >> src >> srcowner >> target >> targetowner >> data;
      atom->add_edge_with_data(src, srcowner, target, targetowner, data);
    } else if (command == 'k') {
      vertex_id_type vid; vertex_color_type color;
      iarc >> vid >> color;
      atom->set_color(vid, color);
    } else if (command == 'l') {
      vertex_id_type vid; uint16_t owner;
      iarc >> vid >> owner;
      atom->set_owner(vid, owner);
    }
  }
  fin.close();
  f.open(filename.c_str(), std::ios::binary | std::ios::app);
}

}
