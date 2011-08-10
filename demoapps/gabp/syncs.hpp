#ifndef _SYNCHS_HPP
#define _SYNCHS_HPP


//helper function to compute the norm between the true solution and the current computed mean
double get_real_norm(const vertex_data& v){
  return pow(v.real - v.cur_mean, 2);
}

//helper function to compute the norm between last round iteration and this round of the current mean
double get_relative_norm(const vertex_data& v){
  return pow(v.cur_mean - v.prev_mean, 2);
}


/**
 * Engine terminates when this returns true
 */
bool termination_condition() {
  double ret = RELATIVE_NORM_KEY.get_val();

  //std::cout<<"I was in term"<<std::endl;
  if (ret < THRESHOLD_KEY.get()){
    std::cout << "Aborting since relative norm is: " << ret << std::endl;
    return true;
  }
  return false;
}


/*
 * aggregate norm of real answer and the current solution
 */
static void apply_func_real(graphlab::any& current_data,
                            const graphlab::any& new_data) {

  double ret = new_data.as<double>();
  ret = sqrt(ret);
  //std::cout << "Real Norm is: " << ret << std::endl;
  // write the final result into the shared data table
  current_data = ret;
}

/*
 * aggregate norm of difference in messages
 */

static void apply_func_relative(graphlab::any& current_data,
                                const graphlab::any& new_data) {

  double ret = new_data.as<double>();
  ret = sqrt(ret);
  std::cout << "Relative Norm is: " << ret << std::endl;
  // write the final result into the shared data table
  current_data = (double)ret;

}




#endif //_SYNCHS_HPP
