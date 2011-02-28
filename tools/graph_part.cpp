#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <fstream>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/command_line_options.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/macros_def.hpp>
using namespace graphlab;
  
typedef distributed_graph<uint32_t, char> graph_type;
double smoothing = 1;

void get_local_count(graph_type &g, 
                      std::vector<size_t> &localpartcounts) {
  foreach(vertex_id_t ownedv, g.owned_vertices()) {
    if (g.get_vertex_data(ownedv) != uint32_t(-1)) {
      localpartcounts[g.get_vertex_data(ownedv)]++;
    }
  }
}
 
bool partition_owned_subgraph(graph_type &g, 
                              std::vector<size_t> &partcounts,
                              std::vector<size_t> &localpartcounts) {
  bool done = true;
  double averagepartcount = 0;
  size_t maxpartsize = 0 ;
  for (size_t i = 0;i < partcounts.size(); ++i) {
    averagepartcount += partcounts[i];
    maxpartsize = std::max<double>(maxpartsize, partcounts[i]);
  }

  averagepartcount = averagepartcount / partcounts.size();
  if (g.rmi.procid() == 0) {
    std::cout << "balance factor: " << maxpartsize << " / " << averagepartcount << " = "
              << maxpartsize / averagepartcount << std::endl;
  }

  averagepartcount += smoothing;
  boost::unordered_set<vertex_id_t> expandedthisround;  
  foreach(vertex_id_t ownedv, g.owned_vertices()) {
    if (g.get_vertex_data(ownedv) != (uint32_t)(-1)) continue;
    // figure out what to join!
    bool recentlyexpanded = false;
    std::map<uint32_t, double> part_to_count;
    // quick check for disconnected vertices
    if (g.in_vertices(ownedv).size() == 0 && g.out_vertices(ownedv).size() == 0) {
      // get the smallest part number
      uint32_t smallestpart = 0;
      uint32_t pcount = uint32_t(-1);
      for (size_t i = 0;i < partcounts.size(); ++i) {
        if (partcounts[i] + localpartcounts[i]  < pcount) { 
          pcount = partcounts[i];
          smallestpart = i;
        }   
      }
      g.set_vertex_data_async(ownedv, smallestpart);
      localpartcounts[smallestpart]++;
      continue; 
    }
    // regular vertices
    foreach(vertex_id_t v, g.in_vertices(ownedv)) {
      if (expandedthisround.find(v) != expandedthisround.end()) {
        recentlyexpanded = true;
        break;
      }
      uint32_t partid = g.get_vertex_data(v);
      if (partid != (uint32_t)(-1)) {
        double partscore = averagepartcount / (partcounts[partid] + localpartcounts[partid] +  smoothing);
        part_to_count[partid] += partscore;
      }
    }
    foreach(vertex_id_t v, g.out_vertices(ownedv)) { 
      if (expandedthisround.find(v) != expandedthisround.end()) {
        recentlyexpanded = true;
        break;
      }
      uint32_t partid = g.get_vertex_data(v);
      if (partid != (uint32_t)(-1)) {
        double partscore = averagepartcount / (partcounts[partid] + localpartcounts[partid] + smoothing);
        part_to_count[partid] += partscore;
      }
    }
    // fast quit
    if (recentlyexpanded || part_to_count.size() == 0) {
      done = false;
      continue;
    }
    
    // find the maximum part
    typedef std::map<uint32_t, double>::value_type ptype;
    uint32_t bestpart = uint32_t(-1);
    double bestpartscore = 0.0;
    foreach(const ptype& p, part_to_count) {
      if (p.second > bestpartscore) {
        bestpartscore = p.second;
        bestpart = p.first;
      }
    }
    if (bestpart == uint32_t(-1)) {
      done = false;
      continue;
    }
    // pick a dd
    expandedthisround.insert(ownedv);
    localpartcounts[bestpart]++;
    g.set_vertex_data_async(ownedv, bestpart);
  }
  return done;
}


bool check_is_done(std::vector<int> vcal) {
  for (size_t i = 0;i < vcal.size(); ++i) {
    if (vcal[i] == false) return false;
  }
  return true;
}

int main(int argc, char** argv) {

  
  global_logger().set_log_level(LOG_DEBUG);
  
  // Parse the command lines
  std::string path;
  std::string outputfile = "parts.txt";
  graphlab::command_line_options clopts("BFS Graph Partitioner", true);
  size_t numparts = 0;
  size_t roundspersync = 1;
  clopts.attach_option("atomindex", &path, path, "atom index file");
  clopts.attach_option("numparts", &numparts, numparts, "numparts");
  clopts.attach_option("outputfile", &outputfile, outputfile, "color output file");
  clopts.attach_option("roundspersync", &roundspersync, roundspersync, "rounds per synchronize");
  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }
  assert(numparts > 0);
  assert(roundspersync >= 1);
  mpi_tools::init(argc, argv);
  dc_init_param param;

  // if not running in DC environment, make atoms
  if (init_param_from_mpi(param) == false) {
    return 0;
  }
  param.initstring="buffered_queued_send=yes";
  distributed_control dc(param);


  graph_type dg(dc, path, true);
  size_t iter = 0;
  dc.barrier();
  // set all vertex data to -1
  foreach(vertex_id_t ownedv, dg.owned_vertices()) {
    dg.vertex_data(ownedv) = uint32_t(-1);
    dg.vertex_is_modified(ownedv);
    dg.increment_vertex_version(ownedv);
  }
  dg.push_all_owned_vertices_to_replicas();
  // synchronize all ghosts
  dg.rmi.full_barrier();
  /// initialization. Select a set of vertices as seeds

  if (dc.procid() == 0) {
    std::set<vertex_id_t> seeds;
    for (size_t i = 0 ; i < numparts; ++i) {
      vertex_id_t v;
      do {
        v = graphlab::random::rand_int(dg.num_vertices() - 1);
      }while(seeds.find(v) != seeds.end());
      seeds.insert(v);
      dg.set_vertex_data(v, i);
    }
  }
  
  dc.full_barrier();
  dc.full_barrier();


  std::vector<size_t> localpartcounts(numparts, 0);
  std::vector<size_t> globalpartcounts(numparts, 1);
  std::vector<std::vector<size_t> > allgathertarget(dc.numprocs()); 
  get_local_count(dg, localpartcounts);
  while(1) {
    iter++;
    bool c = false;
    for (size_t i = 0;i < roundspersync; ++i) {
      c = partition_owned_subgraph(dg, globalpartcounts, localpartcounts);
    }
    
    if (dc.procid() == 0) {
      std::cout << "Check for completion." << std::endl;
    }

    std::vector<int> vcal(dc.numprocs(), 0);
    vcal[dc.procid()] = c;
    dc.all_gather(vcal);

    allgathertarget[dc.procid()] = localpartcounts;
    dc.all_gather(allgathertarget);
    
    std::vector<size_t>  prevglobalpartcounts = globalpartcounts;
    globalpartcounts = std::vector<size_t>(numparts, 0);
    for (size_t i = 0;i <dc.numprocs(); ++i) {
      for (size_t j = 0; j < numparts; ++j) { 
        globalpartcounts[j] += allgathertarget[i][j]; 
      }
    }
    
    bool isdone = check_is_done(vcal);
    if (isdone) break;
    
    
    size_t changedparts = 0;
    for (size_t j = 0; j < numparts; ++j) { 
      changedparts += (globalpartcounts[j] != prevglobalpartcounts[j]);
    }
    if (dc.procid() == 0) {
      for (size_t i = 0; i < std::min<size_t>(numparts, 20); ++i) {
        std::cout << globalpartcounts[i] << " ";
      }
      std::cout << std::endl;
    }
    
    
    
    // if too few partitions changed, we may need to reseed
    if (changedparts < numparts) {
      
      //  partitions are stuck. we need to re-seed a new set of root
      // we seed everyone who has not changed
      size_t meanpart = 0;
      for (size_t i = 0;i < globalpartcounts.size(); ++i) meanpart += globalpartcounts[i];
      meanpart /= dc.numprocs();
      std::vector<vertex_id_t> candidateparts;
      for (size_t i = 0;i < globalpartcounts.size(); ++i) {
        if (globalpartcounts[i] == prevglobalpartcounts[i] && globalpartcounts[i] < meanpart) {
          candidateparts.push_back(i);
        }
      }
      if (dc.procid() == 0) {
        std::cout << "Reseeding " << candidateparts.size() << " partitions" << std::endl;
      }
      if (candidateparts.size()  > 0) {
        std::random_shuffle(candidateparts.begin(), candidateparts.end());
        
        std::vector<std::vector<vertex_id_t> > proposals(dc.numprocs());
        //proposals[dc.procid()].resize(candidateparts.size());
        //every machine proposes a set of vertices for the candidate parts,
        // as well as a count of the number of unassigned vertices
        std::vector<size_t> unassignedvertices(dc.numprocs(), 0);
        foreach(vertex_id_t ownedv, dg.owned_vertices()) {
          if (dg.get_vertex_data(ownedv) == uint32_t(-1)) {
            unassignedvertices[dc.procid()]++;
            if (proposals[dc.procid()].size() < candidateparts.size()) {
              proposals[dc.procid()].push_back(ownedv);
            }
            else {
              size_t sel = random::rand_int(candidateparts.size());
              if (sel < candidateparts.size()) {
                candidateparts[sel] = ownedv;
              }
            }
          }
        }
        dc.all_gather(unassignedvertices, 0);
        // if there are no unassigned vertices then we are done
        size_t allunassigned = 0;
        for (size_t i = 0;i < unassignedvertices.size(); ++i) allunassigned += unassignedvertices[i];
        if (allunassigned > 0) {
          
          dc.gather(proposals, 0);
          if (dc.procid() == 0) {
            //processor 0 now has the task of picking one of the proposed vertices
            for (size_t i = 0;i < candidateparts.size(); ++i) {
              std::vector<vertex_id_t> proposedvertices;
              std::vector<size_t> accweight;
              for (size_t j = 0;j < dc.numprocs(); ++j) {
                if (proposals[j].size() > i) {
                  proposedvertices.push_back(proposals[j][i]);
                  accweight.push_back(unassignedvertices[j]);
                  if (accweight.size() >= 2) accweight[accweight.size() - 1] += accweight[accweight.size() - 2];
                }
              }
              
              if (accweight.size() == 0) continue;
              
              // pick from proposedvertices with weight weight
              size_t r = random::rand_int(accweight[accweight.size() - 1] - 1);
              size_t selv = vertex_id_t(-1);
              for (size_t j = 0; j < accweight.size(); ++j) {
                if (r < accweight[j]) {
                  selv = proposedvertices[j];
                  break;
                }
              }
              assert(selv != vertex_id_t(-1));
              dg.set_vertex_data(selv, candidateparts[i]);
            }
          }
          dc.full_barrier();
          dc.full_barrier(); 
        }
      }
    }
    // check if done
    if (dc.procid() == 0) {
      logstream(LOG_INFO) << "Iteration " << iter << std::endl;
    }
    dg.rmi.full_barrier();
    if (dc.procid() == 0) {
      logstream(LOG_INFO) << "Synchronized " << iter << std::endl;
    }
    dc.barrier();
  }
  
  typedef std::vector<std::pair<vertex_id_t, uint32_t> > vector_of_pairs;
  std::vector<vector_of_pairs> vertex_and_parts(dc.numprocs());
  foreach(vertex_id_t ownedv, dg.owned_vertices()) {
    vertex_and_parts[dc.procid()].push_back(std::make_pair(ownedv, dg.get_vertex_data(ownedv)));
  }
  dc.gather(vertex_and_parts, 0);

  if (dc.procid() == 0) {
    vector_of_pairs result;
    for (size_t i = 0;i < dc.numprocs(); ++i) {
      std::copy(vertex_and_parts[i].begin(), vertex_and_parts[i].end(),
                std::back_inserter(result));
      vector_of_pairs().swap(vertex_and_parts[i]);
      
    }
    std::sort(result.begin(), result.end()); 
    std::ofstream fout;
    fout.open(outputfile.c_str());
    for (size_t i = 0;i < result.size(); ++i) {
      fout << (unsigned int)(result[i].second) << "\n";
    }
  }
  dc.barrier();
  mpi_tools::finalize();
}
