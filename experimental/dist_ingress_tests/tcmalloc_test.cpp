#include <iostream>
#include <vector>
#include <google/malloc_extension.h>


int main(int argc, char** argv) {

   std::vector<int> v;
   for (size_t i = 0; i < 3000*1024; ++i) {
     v.push_back(i);
   }
   size_t hpsize = 0;
   MallocExtension::instance()->GetNumericProperty("generic.heap_size", &hpsize);
   std::cout << "Heap Size: " << hpsize << "\n";

   MallocExtension::instance()->GetNumericProperty("generic.current_allocated_bytes", &hpsize);
   std::cout << "Allocated Size: " << (double)hpsize/(1024*1024) << "MB" << "\n";

} // end of main

