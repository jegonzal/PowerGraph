/*
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

#include <vector>
#include <string>
#include <fstream>

#include <graphlab.hpp>

#include <Python.h>

#define PYFN_TRANSFORMVERTEX 0
#define PYFN_TRANSFORMEDGE   1
#define PYFN_SAVEVERTEX      2
#define PYFN_SAVEEDGE        3
#define PYFN_GATHERAGG       4
#define PYFN_GATHER          5
#define PYFN_APPLY           6
#define PYFN_SCATTER         7
#define PYFN_NEWEDGE         8
#define PYFN_LOADEDGE        9
#define PYFN_STOREEDGE       10
#define PYFN_NEWVERTEX       11
#define PYFN_LOADVERTEX      12
#define PYFN_STOREVERTEX     13
#define PYFN_NEWAGG          14
#define PYFN_LOADAGG         15
#define PYFN_STOREAGG        16

struct {
  const char *fn_name;
  PyObject *fn;
} PyFn[] = {{.fn_name="transformVertex"}, {.fn_name="transformEdge"}, {.fn_name="saveVertex"}, {.fn_name="saveEdge"}, 
           {.fn_name="gatherAgg"}, {.fn_name="gather"}, {.fn_name="apply"}, {.fn_name="scatter"}, 
           {.fn_name="newEdge"}, {.fn_name="loadEdge"}, {.fn_name="storeEdge"}, 
           {.fn_name="newVertex"}, {.fn_name="loadVertex"}, {.fn_name="storeVertex"}, 
           {.fn_name="newAgg"}, {.fn_name="loadAgg"}, {.fn_name="storeAgg"}};

#define PYFN_SIZE  (sizeof(PyFn)/sizeof(PyFn[0]))

PyObject *func_initusermodule;
bool has_edgeclass = false;

#ifndef PYSHARED_LIB
class PythonThreadLocker {
private:  
  PyGILState_STATE state;
public:
  PythonThreadLocker() : state(PyGILState_Ensure()) {}
  ~PythonThreadLocker() { PyGILState_Release(state); }
};
#endif

template<int newmethod_index, int storemethod_index, int loadmethod_index>
class pyobj_class {
  public:
    PyObject *obj;

    pyobj_class() {
#ifndef PYSHARED_LIB
      PythonThreadLocker locker;
#endif
      obj = PyObject_CallObject(PyFn[newmethod_index].fn, NULL);
    }

    pyobj_class(PyObject *no): obj(no) {}

    ~pyobj_class() {
#ifndef PYSHARED_LIB
      PythonThreadLocker locker;
#endif
      Py_XDECREF(obj);
      obj = NULL;
    }

    pyobj_class(const pyobj_class &o) {
#ifndef PYSHARED_LIB
      PythonThreadLocker locker;
#endif
      obj = o.obj;
      Py_XINCREF(obj);
    }

    void operator=(const pyobj_class &o) {
#ifndef PYSHARED_LIB
      PythonThreadLocker locker;
#endif
      Py_XDECREF(obj);
      obj = o.obj;
      Py_XINCREF(obj);
    }

    void save(graphlab::oarchive &oarc) const {
#ifndef PYSHARED_LIB
      PythonThreadLocker locker;
#endif
     
      if (obj == NULL) {
        oarc << std::string("<Py_None>"); 
      } else {
        Py_XINCREF(obj); // we want to keep ownership, prevent SetItem from stealing
        PyObject *pArgs = PyTuple_New(1);
        PyTuple_SetItem(pArgs, 0, obj);       

        PyObject *pValue = PyObject_CallObject(PyFn[storemethod_index].fn, pArgs);
        Py_DECREF(pArgs);
        oarc << std::string(PyString_AsString(pValue));
        Py_DECREF(pValue);
      }     
    }

    void load(graphlab::iarchive &iarc) {
#ifndef PYSHARED_LIB
      PythonThreadLocker locker;
#endif

      std::string s;
      iarc >> s;
      
      if (s == "<Py_None>") {
        Py_XDECREF(obj);
        obj = NULL;
      } else {
        PyObject *pArgs = PyTuple_New(1);
        PyTuple_SetItem(pArgs, 0, PyString_FromString(s.c_str()));
        PyObject *pValue = PyObject_CallObject(PyFn[loadmethod_index].fn, pArgs);
        Py_XDECREF(obj);
        obj = pValue;
        Py_DECREF(pArgs);
      }
    }
};

class agg_class: public pyobj_class<PYFN_NEWAGG, PYFN_STOREAGG, PYFN_LOADAGG> {
  public:
    agg_class(): pyobj_class() {}
    agg_class(PyObject *no): pyobj_class(no) {}

    void operator+=(const agg_class& r) {
#ifndef PYSHARED_LIB
      PythonThreadLocker locker;
#endif

      Py_INCREF(r.obj);  // incref only r.obj and not obj, becuse we'll assign pValue to obj later, prevent SetItem only from stealing r.obj
      PyObject *pArgs = PyTuple_New(2);
      PyTuple_SetItem(pArgs, 0, obj);      
      PyTuple_SetItem(pArgs, 1, r.obj);

      obj = PyObject_CallObject(PyFn[PYFN_GATHERAGG].fn, pArgs);
      Py_DECREF(pArgs);
    }
};

typedef pyobj_class<PYFN_NEWVERTEX, PYFN_STOREVERTEX, PYFN_LOADVERTEX> vertex_data_type;
typedef pyobj_class<PYFN_NEWEDGE, PYFN_STOREEDGE, PYFN_LOADEDGE> edge_data_type;
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;

graph_type *graph;
graphlab::distributed_control *dc;
bool graph_initialized = false;
bool dc_initialized = false;

void transform_vertex(graph_type::vertex_type& vertex) { 
#ifndef PYSHARED_LIB
  PythonThreadLocker locker;
#endif

  PyObject *pArgs = PyTuple_New(1);  // do not incref, obj will be overwritten
  PyTuple_SetItem(pArgs, 0, vertex.data().obj); 
  vertex.data().obj = PyObject_CallObject(PyFn[PYFN_TRANSFORMVERTEX].fn, pArgs);
  Py_DECREF(pArgs);
}

void transform_edge(graph_type::edge_type& edge) {
#ifndef PYSHARED_LIB
  PythonThreadLocker locker;
#endif

  PyObject *pArgs = PyTuple_New(1);  // do not incref, obj will be overwritten
  PyTuple_SetItem(pArgs, 0, edge.data().obj);
  edge.data().obj = PyObject_CallObject(PyFn[PYFN_TRANSFORMEDGE].fn, pArgs);
  Py_DECREF(pArgs);
}

class python_interface:
  public graphlab::ivertex_program<graph_type, agg_class>,
  public graphlab::IS_POD_TYPE {
public:
  agg_class gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
#ifndef PYSHARED_LIB
    PythonThreadLocker locker;
#endif    

    Py_XINCREF(edge.source().data().obj);  // prevent SetItem from stealing references
    Py_XINCREF(vertex.data().obj);
    PyObject *edgedata_arg = (edge.data().obj == NULL) ? Py_None : edge.data().obj;
    Py_XINCREF(edgedata_arg);

    PyObject *pArgs = PyTuple_New(5);
    PyTuple_SetItem(pArgs, 0, edge.source().data().obj);
    PyTuple_SetItem(pArgs, 1, vertex.data().obj);
    PyTuple_SetItem(pArgs, 2, edgedata_arg);
    PyTuple_SetItem(pArgs, 3, PyInt_FromLong(edge.source().num_in_edges()));
    PyTuple_SetItem(pArgs, 4, PyInt_FromLong(edge.source().num_out_edges()));

    PyObject *pValue = PyObject_CallObject(PyFn[PYFN_GATHER].fn, pArgs);
    Py_DECREF(pArgs);
    
    return agg_class(pValue);
  }

  void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
#ifndef PYSHARED_LIB
    PythonThreadLocker locker;
#endif    

    Py_INCREF(total.obj);  // prevent SetItem from stealing total.obj, vertex.obj will be overwriiten
    PyObject *pArgs = PyTuple_New(2);
    PyTuple_SetItem(pArgs, 0, vertex.data().obj);
    PyTuple_SetItem(pArgs, 1, total.obj);

    vertex.data().obj = PyObject_CallObject(PyFn[PYFN_APPLY].fn, pArgs);
    Py_DECREF(pArgs);
  }

  void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
#ifndef PYSHARED_LIB
    PythonThreadLocker locker;
#endif    

    Py_XINCREF(edge.source().data().obj);  // prevent SetItem from stealing references
    Py_XINCREF(vertex.data().obj);
    PyObject *edgedata_arg = (edge.data().obj == NULL) ? Py_None : edge.data().obj;
    Py_XINCREF(edgedata_arg);


    PyObject *pArgs = PyTuple_New(5);
    PyTuple_SetItem(pArgs, 0, edge.source().data().obj);
    PyTuple_SetItem(pArgs, 1, vertex.data().obj);
    PyTuple_SetItem(pArgs, 2, edgedata_arg);
    PyTuple_SetItem(pArgs, 3, PyInt_FromLong(edge.source().num_in_edges()));
    PyTuple_SetItem(pArgs, 4, PyInt_FromLong(edge.source().num_out_edges()));

    PyObject *pValue = PyObject_CallObject(PyFn[PYFN_SCATTER].fn, pArgs);
    Py_DECREF(pArgs);
    
    int signal_result = PyInt_AsLong(PyTuple_GetItem(pValue, 0));
    if (signal_result) {
      context.signal(edge.target());
    }
      
    PyObject *edgedata_result = PyTuple_GetItem(pValue, 1);
    if (edgedata_result != Py_None) {
      Py_DECREF(edge.data().obj);
      edge.data().obj = edgedata_result;
    }
    
    PyObject *vertexdata_result = PyTuple_GetItem(pValue, 2);
    if (vertexdata_result != Py_None) {
      Py_DECREF(edge.target().data().obj);
      edge.target().data().obj = vertexdata_result;
    }

    Py_DECREF(pValue);
  }
};

struct py_writer {
  
  std::string save_vertex(graph_type::vertex_type v) {
#ifndef PYSHARED_LIB
    PythonThreadLocker locker;
#endif

    Py_INCREF(v.data().obj);
    PyObject *pArgs = PyTuple_New(1);
    PyTuple_SetItem(pArgs, 0, v.data().obj);
    PyObject *pValue = PyObject_CallObject(PyFn[PYFN_SAVEVERTEX].fn, pArgs);
    Py_DECREF(pArgs);

    std::stringstream strm;
    strm << v.id() << "\t" << PyString_AsString(pValue) << "\n";
    Py_DECREF(pValue);
    return strm.str();
  }
  
  std::string save_edge(graph_type::edge_type e) { 
#ifndef PYSHARED_LIB
    PythonThreadLocker locker;
#endif

    Py_INCREF(e.data().obj);
    PyObject *pArgs = PyTuple_New(1);
    PyTuple_SetItem(pArgs, 0, e.data().obj);
    PyObject *pValue = PyObject_CallObject(PyFn[PYFN_SAVEEDGE].fn, pArgs);
    Py_DECREF(pArgs);

    std::stringstream strm;
    strm << PyString_AsString(pValue) << "\n";
    Py_DECREF(pValue);
    return strm.str();
  }
};

int init_python(const char *python_script) {
#ifndef PYSHARED_LIB
  // Initialize Python
  Py_Initialize();
  PyEval_InitThreads();

  PyThreadState* tcur = PyThreadState_Get();
  PyThreadState_Swap(NULL);
  PyThreadState_Clear(tcur);
  PyThreadState_Delete(tcur);

  PyEval_ReleaseLock();
#endif

{
#ifndef PYSHARED_LIB
  PythonThreadLocker locker;
#endif

  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append('')");

  PyObject *pNameWrap = PyString_FromString("wrappers");
  PyObject *pModuleWrap = PyImport_Import(pNameWrap);
  Py_DECREF(pNameWrap);

  if (pModuleWrap == Py_None) {
    PyErr_Print();
    dc->cout() << "Failed to load wrapper script\n";
    return EXIT_FAILURE;
  }

  func_initusermodule = PyObject_GetAttrString(pModuleWrap, "initUserModule");

  for (int i = 0; i < PYFN_SIZE; i++) {
    PyFn[i].fn = PyObject_GetAttrString(pModuleWrap, PyFn[i].fn_name);
    if (PyFn[i].fn == Py_None || PyFn[i].fn == NULL) {
      dc->cout() << "Cannot load python function " << PyFn[i].fn_name << " from wrapper module\n";
      return EXIT_FAILURE;
    }
  }
  
  Py_DECREF(pModuleWrap);

  PyObject *pArgs = PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, PyString_FromString(python_script));
  PyObject *pValue = PyObject_CallObject(func_initusermodule, pArgs);
  Py_DECREF(pArgs);
  has_edgeclass = pValue == Py_True;
  Py_DECREF(pValue);
}
  
  return EXIT_SUCCESS;
}

void init() {
  if (dc_initialized) {
    dc->cout() << "DC already initialized\n";
    return;
  }

  int argc = 0;
  char **argv = NULL;
  
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  dc = new graphlab::distributed_control();
  global_logger().set_log_level(LOG_INFO);
  dc_initialized = true;
}

int init_graph(const char *graph_dir, const char *format, const char *python_script) {
  if (!dc_initialized) {                      
    dc->cout() << "DC not initialized\n";
    return EXIT_FAILURE;
  }
  if (graph_initialized) {
    dc->cout() << "Graph already initialized\n";
    return EXIT_SUCCESS;
  }

  if (graph_dir == NULL) {
    dc->cout() << "Graph not specified. Cannot continue";
    return EXIT_FAILURE;
  }

  if (init_python(python_script)) {
    dc->cout() << "Python initialization failure.";
    return EXIT_FAILURE;
  }

  graphlab::command_line_options clopts("Python algorithm.");
  clopts.set_ncpus(1);

  // Build the graph ----------------------------------------------------------
  graph = new graph_type(*dc, clopts);
  dc->cout() << "Loading graph in format: "<< format << std::endl;
  graph->load_format(std::string(graph_dir), std::string(format));
  // must call finalize before querying the graph
  graph->finalize();
  dc->cout() << "#vertices: " << graph->num_vertices() << " #edges: " << graph->num_edges() << std::endl;
  graph_initialized = true;

  return EXIT_SUCCESS;
}

int transform_graph() {
  if (!dc_initialized) {
    dc->cout() << "DC not initialized\n";
    return EXIT_FAILURE;
  }
  if (!graph_initialized) {
    dc->cout() << "Graph not initialized\n";
    return EXIT_FAILURE;
  }

//  graph->transform_edges(transform_edge);
  graph->transform_vertices(transform_vertex);
  return EXIT_SUCCESS;
}

int gas_graph(const char *exec_type) {
  if (!dc_initialized) {                      
    dc->cout() << "DC not initialized\n";
    return EXIT_FAILURE;
  }
  if (!graph_initialized) {
    dc->cout() << "Graph not initialized\n";
    return EXIT_FAILURE;
  }

  // Running The Engine -------------------------------------------------------
  graphlab::command_line_options clopts("Python algorithm.");
  clopts.set_ncpus(1);

  graphlab::omni_engine<python_interface> engine(*dc, *graph, std::string(exec_type), clopts);
  engine.signal_all();
  engine.start();
  const float runtime = engine.elapsed_seconds();
  dc->cout() << "Finished Running engine in " << runtime << " seconds." << std::endl;
  return EXIT_SUCCESS;
}

int save_graph(const char *save_prefix, const int use_gzip, const int save_vertices, const int save_edges) {
  if (!dc_initialized) {                      
    dc->cout() << "DC not initialized\n";
    return EXIT_FAILURE;
  }
  if (!graph_initialized) {
    dc->cout() << "Graph not initialized\n";
    return EXIT_FAILURE;
  }

  graph->save(std::string(save_prefix), py_writer(), use_gzip, save_vertices, save_edges);  
  return EXIT_SUCCESS;
}

int done_graph() {
  if (!dc_initialized) {                      
    dc->cout() << "DC not initialized\n";
    return EXIT_FAILURE;
  }
  if (!graph_initialized) {
    dc->cout() << "Graph already done\n";
    return EXIT_SUCCESS;
  }
  delete graph;
  graph_initialized = false;

  return EXIT_SUCCESS;
}

void done() {
  if (!dc_initialized) {                      
    dc->cout() << "DC already done\n";
    return;
  }
  if (graph_initialized) {
    dc->cout() << "Graph still initialized\n";
    return;
  }
  graphlab::mpi_tools::finalize();
  delete dc;
  dc_initialized = false;
}

int main(int argc, char **argv) {
  graphlab::command_line_options clopts("Python algorithm.");

  std::string python_script = "";
  std::string graph_dir = "";
  std::string format = "";
  std::string exec_type = "synchronous";
  std::string save_prefix = "";

  clopts.attach_option("script", python_script, "Python script. Required ");
  clopts.attach_option("graph", graph_dir, "The graph file. Required ");
  clopts.add_positional("graph");
  clopts.attach_option("format", format, "The graph file format");
  clopts.attach_option("engine", exec_type, "The engine type synchronous or asynchronous");
  clopts.attach_option("saveprefix", save_prefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");

  if(!clopts.parse(argc, argv)) {
    dc->cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  init();
  init_graph(graph_dir.c_str(), format.c_str(), python_script.c_str());
  transform_graph();
  gas_graph(exec_type.c_str());
  save_graph(save_prefix.c_str(), 0, 1, 0);
  done_graph();
  done();

  return EXIT_SUCCESS;
} 



