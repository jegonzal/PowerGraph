%module py_graphlab

%{
void init();
int init_graph(const char *graph_dir, const char *format, const char *python_script);
int transform_graph();
int gas_graph(const char *exec_type);
int save_graph(const char *save_prefix, const int use_gzip, const int save_vertices, const int save_edges);
int done_graph();
void done();
%}

%init %{
    init();
    atexit(done);
%}

void init();
int init_graph(const char *graph_dir, const char *format, const char *python_script);
int transform_graph();
int gas_graph(const char *exec_type);
int save_graph(const char *save_prefix, const int use_gzip, const int save_vertices, const int save_edges);
int done_graph();
void done();

