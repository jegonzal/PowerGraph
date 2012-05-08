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
 
/**
 * @file message.cpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
 
#include "json_message.hpp"
#include <graphlab/macros_def.hpp>

using namespace graphlab;
namespace json = rapidjson;

///////////////////////////////// C CALLBACKS //////////////////////////////////

/**
 * If parent is an object, adds child to pending key and removes pending key;
 * If parent is an array, adds child to end of array
 * Otherwise, throws error.
 */
static json::Value &
add_to_parent(json::Value& parent,
              json::Value& child,
              json::Document::AllocatorType& allocator,
              struct jsonsl_state_st *state){

  if (parent.IsArray()) {
      json::Value &val = parent.PushBack(child, allocator);
      state->data = &val;
      return val;
  }
  
  if (parent.IsObject()) {
    json::Value &val = parent.AddMember(parent["pending-key"], child, allocator);
    state->data = &val;
    parent.RemoveMember("pending-key");
    return val;
  }

  logstream(LOG_ERROR) << parent.GetString() << std::endl;
  logstream(LOG_ERROR) << "Requested to add to non-container parent type!" << std::endl;
  CHECK(false); // TODO error handling
    
}

/**
 * Called to create new JSON element
 */
static void
create_new_element(jsonsl_t jsn,
                   jsonsl_action_t action,
                   struct jsonsl_state_st *state,
                   const char *buf){

    logstream(LOG_DEBUG) << "pushed with " << buf << std::endl;

    json::Value child;
    json::Value *parent_ptr = NULL;
    
    // extract state from stack
    CHECK(state);
    struct jsonsl_state_st *last_state = jsonsl_last_state(jsn, state);
    if (NULL == last_state){
      // root element
      CHECK(JSONSL_T_OBJECT == state->type);
      state->data = jsn->data;
      return;
    }

    parent_ptr = (json::Value *)last_state->data;
    json::Document *document = (json::Document *) jsn->data;
    json::Document::AllocatorType& allocator = document->GetAllocator();
    
    CHECK(parent_ptr);
    json::Value &parent = *parent_ptr;

    switch(state->type) {
    case JSONSL_T_SPECIAL: // TODO!
    case JSONSL_T_STRING: {
        return; // wait for complete string
    }
    case JSONSL_T_HKEY: {
        return; // wait for complete string
    }
    case JSONSL_T_LIST:
        child.SetArray();
        break;
    case JSONSL_T_OBJECT:
        child.SetObject();
        break;
    default:
        logstream(LOG_ERROR) << "Unhandled type: " << state->type << std::endl;
        CHECK(false);
        // TODO error handling
        break;
    }

    add_to_parent(parent, child, allocator, state);
    
}

static void
end_element (jsonsl_t jsn,
            jsonsl_action_t action,
            struct jsonsl_state_st *state,
            const char *buf){
            
  logstream(LOG_DEBUG) << "popped at " << buf << std::endl;
  CHECK(state);
  
  struct jsonsl_state_st *last_state = jsonsl_last_state(jsn, state);
  if (NULL == last_state){
    // root element
    jsn->data = NULL;
    return;
  }

  json::Value *parent_ptr = (json::Value *)last_state->data;
  json::Document *document = (json::Document *) jsn->data;
  json::Document::AllocatorType& allocator = document->GetAllocator();
  
  CHECK(parent_ptr);
  json::Value &parent = *parent_ptr;
  
  // complete string (add null, copy)
  if ('"' == *buf){
  
    *(char *)buf = '\0';
    
    switch(state->type){
      case JSONSL_T_HKEY: {
        json::Value child(buf-state->pos_begin+1, allocator);
        json::Value &val = parent.AddMember("pending-key", child, allocator);
        // not copied: this is a bug - lexer requires complete strings
        state->data = &val;
        break;
      }
      case JSONSL_T_STRING: {
        json::Value child(buf-state->pos_begin+1, allocator);
        add_to_parent(parent, child, allocator, state);
        break;
      }
      default: break;   
    }
    
  }
  
}

/**
 * Called on JSON parse error.
 */
static int error_callback(jsonsl_t jsn,
                jsonsl_error_t err,
                struct jsonsl_state_st *state,
                char *at){
    logstream(LOG_ERROR) << "Got error at pos " << jsn->pos << ": " << jsonsl_strerror(err) << std::endl;
    if (at) logstream(LOG_ERROR) << "Remaining text: " << std::string(at) << std::endl;
    // TODO error handling
    return 0;
}

/////////////////////////////// INSTANCE MEMBERS ///////////////////////////////

json_message::
  json_message(const std::string method, const std::string state) :
  mdocument(), mjsn(jsonsl_new(MAX_LEVELS)) {

  // make object root
  mdocument.SetObject();
  
  // add method if exists
  if (method.length() > 0){
  	json::Value methodv(method.c_str(), mdocument.GetAllocator());
  	mdocument.AddMember("method", methodv, mdocument.GetAllocator());
  }
  
  // add state if exists
  if (state.length() > 0){
  	json::Value statev (state.c_str(), mdocument.GetAllocator());
	  mdocument.AddMember("state", statev, mdocument.GetAllocator());
	}
	
	// TODO: optimize jsn -> how many levels do we need?
	jsonsl_enable_all_callbacks(mjsn);
	mjsn->data = &mdocument;
  mjsn->action_callback = NULL;
  mjsn->action_callback_PUSH = create_new_element;
  mjsn->action_callback_POP = end_element;
  mjsn->error_callback = error_callback;
  mjsn->max_callback_level = 20;
  
}

bool json_message::feed(byte *data, std::size_t nread){
  jsonsl_feed(mjsn, data, nread);
  return (NULL != mjsn->data);
}

json_message::
  ~json_message(){
  jsonsl_destroy(mjsn);
}

///////////////////////////////// CLASS MEMBERS ////////////////////////////////
const int json_message::MAX_LEVELS = 0x1000;

#include <graphlab/macros_undef.hpp>
