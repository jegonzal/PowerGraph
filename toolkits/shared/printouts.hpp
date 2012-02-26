#ifndef PRINTOUTS
#define PRINTOUTS

void print_vec(const char * name, DistVec & vec, bool high = false){
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
void print_vec(const char * name, vec & pvec, bool high = false){
 int i;
 printf("%s[%d]\n", name, 0);
 for (i= 0; i< pvec.size(); i++){
  if (high)
   printf("%15.15lg ", pvec[i]);
  else
   printf("%.5lg ", pvec[i]);
  }
 printf("\n");
}


#define PRINT_VEC(a) print_vec(#a,a,0)
#define PRINT_VEC2_HIGH(a,i) print_vec(#a,a[i],1)
#define PRINT_INT(a) printf("%s: %d\n", #a, a);
#define PRINT_DBL(a) printf("%s: %.5lg\n", #a, a);

#endif
