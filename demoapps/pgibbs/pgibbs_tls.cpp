/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#include "pgibbs_tls.hpp"

pthread_key_t pgibbs_tls_key;

pgibbs_tls* create_pgibbs_tls() {
  ASSERT_EQ(pthread_getspecific(pgibbs_tls_key), NULL);
  pgibbs_tls* data = new pgibbs_tls();
  ASSERT_NE(data, NULL);
  pthread_setspecific(pgibbs_tls_key, data);
  return data;
}

pgibbs_tls& get_pgibbs_tls() {
  pgibbs_tls* tls =
    reinterpret_cast<pgibbs_tls*>
    (pthread_getspecific(pgibbs_tls_key) );
  // If no tsd be has been associated, create one
  if(tls == NULL) tls = create_pgibbs_tls();
  ASSERT_NE(tls, NULL);
  return *tls;
}

void destroy_pgibbs_tls(void* ptr) {
  pgibbs_tls* tls = 
    reinterpret_cast<pgibbs_tls*>(ptr);
  if(tls != NULL) delete tls;

}


struct pgibbs_tls_key_creater {
  pgibbs_tls_key_creater( )  {
    pthread_key_create(&pgibbs_tls_key,
                       destroy_pgibbs_tls);
  }
};
static const pgibbs_tls_key_creater make_pgibbs_tls_key;


