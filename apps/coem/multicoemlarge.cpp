/*
 * \author akyrola
 *  Starter for CoEM
 */

#include <cstdio>
#include <string>

#include <graphlab/logger/logger.hpp>

#include "multicoemapp.hpp"


int main(int argc,  const char *argv[]) {
  global_logger().set_log_level(0);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "MultiCoEM / ReadTheWeb starting\n");
  
  
  multicoemapp * app = new multicoemapp(
					"/mnt/bigbrofs/usr5/graphlab/testdata/coem/justin/large/cat_nps",
					"/mnt/bigbrofs/usr5/graphlab/testdata/coem/justin/large/cat_contexts", 
					"/mnt/bigbrofs/usr5/graphlab/testdata/coem/justin/large/aapo-ctx-idx.txt", 
					"/mnt/bigbrofs/usr5/graphlab/testdata/coem/justin/large/seeds/", 
					"/mnt/bigbrofs/usr5/graphlab/testdata/coem/justin/large/seeds-neg/");
  app->parse_args(argc, (char**) argv);
  app->prune_edges = true;
  app->start();
}
