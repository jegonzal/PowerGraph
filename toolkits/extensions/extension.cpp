#include <map>
#include <string>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/lockfree_push_back.hpp>
#include "extension_gas.hpp"

namespace graphlab {
namespace extension {

var& operator+=(var& value, const var& other) {
  const double* other_double = boost::get<double>(&other);
  const std::string* other_string = boost::get<std::string>(&other);
  const Eigen::VectorXd* other_vector = boost::get<Eigen::VectorXd>(&other);
  const Eigen::MatrixXd* other_matrix = boost::get<Eigen::MatrixXd>(&other);

  if ( double* val = boost::get<double>( &value) ) {
    if (other_double!= NULL) value = (double)(*val) + (*other_double);
    else ASSERT_MSG(false, "Type mismatch in operator+=");
  } else if ( std::string* val = boost::get<std::string>( &value) ) {
    if (other_string != NULL) (*val) += (*other_string);
    else ASSERT_MSG(false, "Type mismatch in operator+=");
  } else if ( Eigen::VectorXd* val = boost::get<Eigen::VectorXd>( &value) ) {
    if (other_vector!= NULL) (*val) += (*other_vector);
    else ASSERT_MSG(false, "Type mismatch in operator+=");
  } else if ( Eigen::MatrixXd* val = boost::get<Eigen::MatrixXd>( &value) ) {
    if (other_matrix!= NULL) (*val) += (*other_matrix);
    else ASSERT_MSG(false, "Type mismatch in operator+=");
  }
  return value;
}









// lets get more than we will ever need so it will never need to resize
std::vector<gas_op_descriptor> descriptor_set(65536);
lockfree_push_back<std::vector<gas_op_descriptor> > 
          descriptor_access(descriptor_set, 0);
 


var vars::empty_var;
}
}
