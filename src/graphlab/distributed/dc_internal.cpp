#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/distributed/dc_internal.hpp>
namespace graphlab {

void destroy_dc_thread_local_struct(void* ptr) ;

struct dc_thread_keys {
  pthread_key_t GRAPHLAB_DC_BUFFER_ID;
  dc_thread_keys() : GRAPHLAB_DC_BUFFER_ID(1){ 
    pthread_key_create(&GRAPHLAB_DC_BUFFER_ID,
                        destroy_dc_thread_local_struct);
  }
};

static const dc_thread_keys dckeys;



dc_thread_local_struct* create_thread_dc_buffer(size_t thread_id = 0) {
  // Require that the data not yet exist
  assert(pthread_getspecific(dckeys.GRAPHLAB_DC_BUFFER_ID) == NULL);
  // Create the data
  dc_thread_local_struct* data = new dc_thread_local_struct;
  ASSERT_NE(data, NULL);
  memset(data, 0, sizeof(dc_thread_local_struct));
  // Set the data
  pthread_setspecific(dckeys.GRAPHLAB_DC_BUFFER_ID, data);
  // Return the associated tsd
  return data;
}


dc_thread_local_struct& get_thread_dc_buffer() {
  // get the tsd
  dc_thread_local_struct* tsd =
    reinterpret_cast<dc_thread_local_struct*>
    (pthread_getspecific(dckeys.GRAPHLAB_DC_BUFFER_ID));
  // If no tsd be has been associated, create one
  if(tsd == NULL) tsd = create_thread_dc_buffer();
  assert(tsd != NULL);
  return *tsd;
}


void destroy_dc_thread_local_struct(void* ptr) {
  dc_thread_local_struct* tsd =
      reinterpret_cast<dc_thread_local_struct*>
      (pthread_getspecific(dckeys.GRAPHLAB_DC_BUFFER_ID));

  if (tsd != NULL && tsd->sendbuffer != NULL) free(tsd->sendbuffer);
}

}