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
  

boost::posix_time::ptime myEpoch(boost::gregorian::date(1970,boost::gregorian::Jan,1)); // Or whatever your epocj is.

unsigned long int datestr2uint64(const std::string & data){

//std::string ts("2002-01-20 23:59:59.000");
//ptime t(time_from_string(ts))
  std::string ts = "20" + data.substr(4,2) + "-" + data.substr(2,2) + "-" + data.substr(0,2) +
                   " " + data.substr(7,2) + ":" + data.substr(9,2) + ":" + data.substr(11,2) + ".000";
  boost::posix_time::ptime t(boost::posix_time::time_from_string(ts));
  boost::posix_time::time_duration myTimeFromEpoch = t - myEpoch;
  return myTimeFromEpoch.ticks();
}
