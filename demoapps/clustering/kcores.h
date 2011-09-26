#ifndef _KCORES_H
#define _KCORES_H


struct kcores_data{
  bool active;
  int kcore;
  int degree;

  kcores_data(){
    active = true;
    kcore = -1;
    degree = 0;
  }

   
  void save(graphlab::oarchive& archive) const{
    archive << active << kcore << degree;
  } 
   
  void load(graphlab::iarchive& archive){
    archive >> active >> kcore >> degree;
  } 


};

struct kcores_edge{
  float weight;

  kcores_edge(){
    weight = 1;
  };

 void save(graphlab::oarchive& archive) const{
   archive << weight;
 } 
   
  void load(graphlab::iarchive& archive){
   archive >> weight;
  } 

};




#endif //_KCORES_H
