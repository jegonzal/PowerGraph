#ifndef PRINTOUTS
#define PRINTOUTS

void print_vec(const char * name, const DistVec & vec, bool high = false){
 int i;
 printf("%s[%d]\n", name, vec.offset);
 for (i=vec.start; i< vec.end; i++){
  if (high)
   printf("%15.15lg ", pgraph->vertex_data(i).pvec[vec.offset]);
  else
   printf("%.5lg ", pgraph->vertex_data(i).pvec[vec.offset]);
  }
 printf("\n");
}
void print_vec(const char * name, const vec & pvec, bool high = false){
 printf("%s\n", name);
 for (int i= 0; i< pvec.size(); i++){
  if (high)
   printf("%15.15lg ", pvec[i]);
  else
   printf("%.5lg ", pvec[i]);
  }
 printf("\n");
}
void print_mat(const char * name, const mat & pmat, bool high = false){
 printf("%s\n", name);
 for (int i= 0; i< pmat.rows(); i++){
  for (int j=0; j< pmat.cols(); j++){
    if (high)
      printf("%15.15lg ", get_val(pmat, i, j));
    else
     printf("%.5lg ", get_val(pmat, i, j));
  }
  printf("\n");
  
  }
}

void print_vec_pos(std::string name, vec & v, int i){
  if (i == -1)
    printf("%s\n", name.c_str());
  else {
    printf("%s[%d]: %.5lg\n", name.c_str(), i, v[i]);
    return;
  }
  for (int j=0; j< v.size(); j++){
   printf("%.5lg", v(j));
   if (v.size() > 1)
    printf(" ");
  }
  printf("\n");
}


#define PRINT_VEC(a) print_vec(#a,a,0)
#define PRINT_VEC2(a,b) print_vec(a,b,0)
#define PRINT_VEC3(a,b,c) print_vec_pos(a,b,c)
#define PRINT_VEC2_HIGH(a,i) print_vec(#a,a[i],1)
#define PRINT_INT(a) printf("%s: %d\n", #a, a);
#define PRINT_NAMED_INT(a,b) printf("%s: %d\n",a, b);
#define PRINT_DBL(a) printf("%s: %.5lg\n", #a, a);
#define PRINT_NAMED_DBL(a,b) printf("%s: %.5lg\n", a, b);
#define PRINT_MAT(a) print_mat(#a, a, 0);
#define PRINT_MAT2(a,b) print_mat(a,b,0);
#endif
