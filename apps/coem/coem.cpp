/*
 * \author akyrola
 *  Starter for CoEM
 */

#include <cstdio>
#include <logger/logger.hpp>

#include "coemapp.hpp"


int main(int argc,  const char *argv[]) {
  global_logger().set_log_level(0);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "CoEM / ReadTheWeb starting\n");
  
  // Defaults
  if (argc == 1) {
    argc = 6;   
    argv = (const char**) malloc(sizeof(char*)*7);
    int k = 1;
    // THIS IS GENERATES COMPILER WARNINGS:
    // FIX:
    argv[k++] = "../testdata/coem/nps.txt";
    argv[k++] = "../testdata/coem/contexts.txt";
    argv[k++] = "../testdata/coem/matrix.txt";
    argv[k++] = "../testdata/coem/city-seeds.txt";
    argv[k++] = "../testdata/coem/city-neg-seeds.txt";
    
    logger(LOG_INFO, "*** Using default files");
  }
  
  if (argc < 6) {
    logger(LOG_ERROR, 
           "Usage: coem [nps-file] [context-file] "
           "[matrix-file] [seeds-file] [negseeds-file] "
           "[patopms]\n");
    return 1;
  }
  
  coemapp * app = new coemapp(argv[1], argv[2], argv[3], argv[4], argv[5]);
  app->parse_args(argc, (char**) argv);
  app->start();
}
