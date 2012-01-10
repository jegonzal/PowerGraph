/*
YXVaVQJfYZp BqFnHyiRwam 050803 235959 28
xtGBWvpYgYK jdJBsGbXUwu 050803 235959 242
ZpYQLFFKyTa atslZokWZRL 050803 235959 504
WMjbYygLglR BqFnCfuNgio 050803 235959 51
hcLiEoumskU RcNSJSEmidT 050803 235959 7
qBvUQlMABPv atslBvPNusB 050803 235959 3609
jdSqVjxPlBn BqFnHyiRwam 050803 235959 23
VCWyivVLzRr atslSbOOWXz 050803 235959 8411
PnLaFqLJrEV atslZokWZRL 050803 235959 8806
PnLaCsEnqei atslBvPNusB 050803 235959 590
*/
#include "boost/date_time/posix_time/posix_time.hpp"
  

boost::posix_time::ptime myEpoch(boost::gregorian::date(2005,boost::gregorian::Aug,1)); // Or whatever your epocj is.

struct cached_time{
std::string cached_data;
unsigned long int cached_timeret;
unsigned long int cached_dateret;
};

cached_time cache[16];

//input string in the format: 050803 235959
unsigned long int datestr2uint64(const std::string & data, int & timeret, int & dateret, int thread_id){

  if (cache[thread_id].cached_data == data){
    timeret = cache[thread_id].cached_timeret;
    dateret = cache[thread_id].cached_dateret;
    return -1;
  } 

//std::string ts("2002-01-20 23:59:59.000");
//ptime t(time_from_string(ts))
  std::string ts = "20" + data.substr(0,2) + "-" + data.substr(2,2) + "-" + data.substr(4,2) +
                   " " + data.substr(7,2) + ":" + data.substr(9,2) + ":" + data.substr(11,2) + ".000";
  boost::posix_time::ptime t(boost::posix_time::time_from_string(ts));
  boost::posix_time::time_duration myTimeFromEpoch = t - myEpoch;
  unsigned long int ticks = myTimeFromEpoch.ticks() / 1000000;
  timeret = ticks % (3600*24);
  dateret = ticks / (3600*24);
  cache[thread_id].cached_timeret = timeret;
  cache[thread_id].cached_dateret = dateret;
  cache[thread_id].cached_data = data;
  return myTimeFromEpoch.ticks();
}
