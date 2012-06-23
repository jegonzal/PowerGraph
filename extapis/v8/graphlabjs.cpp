#define OUTPUTLEVEL   LOG_INFO

#include <graphlab.hpp>
#include <v8.h>
#include "cvv8.hpp"

using namespace v8;
using namespace graphlab;
namespace cv = cvv8;

namespace graphlab {

  // TODO: refactor graph operations to another graph object
  class pilot {
  private:
    
    // TODO: allow more generic types
    typedef double vertex_data_type;
    typedef double edge_data_type;
    typedef distributed_graph<vertex_data_type, edge_data_type> graph_type;

    static graphlab_options opts;
    distributed_control dc;
    graph_type graph;

  public: 
    
    // TODO: how many distributed controls can a pilot have in his lifetime? 
    pilot() : dc(), graph(dc, opts) {}
    
    void ping(){ std::cout << "pong." << std::endl; }
 
    // TODO: how many graphs can a pilot have in his lifetime? 
    void load_graph(const std::string &path, const std::string &format){
      graph.load_format(path, format);
      graph.finalize();
    }

    void fly(const Handle<Function> &function){
      function->Call(function, 0, NULL);
    }

    /**
     * Adds a JS binding of the class to the given object. Throws
     * a native exception on error.
     */
    static void setup_bindings(v8::Handle<v8::Object> const &dest){
      cv::ClassCreator<pilot>::Instance().SetupBindings(dest);
    }

    /**
     * Saves command line options for this session.
     */
    static void set_clopts(const graphlab_options &clopts){ opts = clopts; }

  };

};

graphlab_options pilot::opts;

namespace cvv8 {

    template<>
    CVV8_TypeName_IMPL((pilot), "pilot");

    // A helper to support converting from pilot to its JS handle.
    typedef NativeToJSMap<pilot> BMap;

    // pilot Ctors we want to bind to v8 (there are several other ways to do this):
    typedef Signature<pilot (
      CtorForwarder<pilot *()>
    )> pilotCtors;

   /**
      This policy class is required unless you just want to bind to 
      the default constructor. It creates native objects for the 
      underlying binding code.
    */
    template <>
    class ClassCreator_Factory<pilot> {
    public:
        typedef pilot * ReturnType;
        static ReturnType Create( v8::Persistent<v8::Object> & jsSelf, v8::Arguments const & argv ){
          typedef cv::CtorArityDispatcher<pilotCtors> Proxy;
          pilot * b = Proxy::Call( argv );
          if( b ) BMap::Insert( jsSelf, b );
          return b;
        }
        static void Delete( ReturnType obj ){
          BMap::Remove( obj );
          delete obj;
        }
    };

    /** Disable subclassing from within javascript */
    template <>
    struct ClassCreator_SearchPrototypeForThis<pilot> : Opt_Bool<false> {};
    
    /**
       Required specialization so that the conversion API can derive
       the native 'this' object from v8::Arguments::This() and from
       function arguments of our bound type.

       This implementation works by using the plumbing installed
       by ClassCreator.
    */
    template <>
    struct JSToNative<pilot>
        : JSToNative_ClassCreator<pilot> {};

    /**
        Native-to-JS conversion. This conversion is only possible when
        we explicitly add support to the class-binding code to add
        the necessary binding metadata. Alternately, if the bound class
        contains v8-defined data types, e.g. a Handle<Object> referring
        to itself then implementing NativeToJS is easy to do - just
        return the handle held by the native.
        
        In this example we're using the NativeToJSMap helper code to 
        plug in/unplug our bindings during native object 
        construction/destruction via the the ClassCreator's factory 
        policy (ClassCreator_Factory<pilot>).        
    */
    template <>
    struct NativeToJS<pilot>
        :  NativeToJSMap<pilot>::NativeToJSImpl {};
 

    template <>
    struct ClassCreator_SetupBindings<pilot> {

      static void Initialize(Handle<Object> const & dest) {

        logstream(LOG_INFO) << "== Preparing cockpit " << std::flush;

        ////////////////////////////////////////////////////////////
        // Bootstrap class-wrapping code...
        typedef ClassCreator<pilot> CC;
        CC & cc( CC::Instance() );
        if( cc.IsSealed() ) {
          cc.AddClassTo( TypeName<pilot>::Value, dest );
          logstream(LOG_INFO) << "== cockpit ready ==>" << std::endl; 
          return;
        }

        ////////////////////////////////////////////////////////////
        // Bind some member functions...
        logstream(LOG_INFO) << "== strapping pilot " << std::flush;
        cc("ping", MethodToInCa<pilot, void (), &pilot::ping>::Call)
          ("destroy", CC::DestroyObjectCallback)
          ("loadGraph", 
            MethodToInCa<pilot, void (const std::string&, const std::string&),
              &pilot::load_graph>::Call)
          ("fly",
            MethodToInCa<pilot, void (const Handle<Function> &),
              &pilot::fly>::Call);

        ////////////////////////////////////////////////////////////
        // Add class to the destination object...
        cc.AddClassTo( TypeName<pilot>::Value, dest );

        // sanity checking. This code should crash if the basic stuff is horribly wrong
        logstream(LOG_INFO) << "== running test flight " << std::flush;
        Handle<Value> vinst = cc.NewInstance (0,NULL);
        pilot * native = CastFromJS<pilot>(vinst);
        if (0 == native)
          logstream(LOG_FATAL) << "== BROKEN. Unable to instantiate test flight." << std::endl;

        logstream(LOG_INFO) << "== success ==>"
          << std::endl;
        CC::DestroyObject( vinst );

        logstream(LOG_EMPH) << "Cockpit ready." << std::endl;

     }

  };

};

/**
 * Initializes v8 shell and setups bindings.
 */
static int v8_main(const std::string& script){

  cv::Shell shell;
  shell.SetupDefaultBindings();
  HandleScope scope;
  
  pilot::setup_bindings(shell.Global());
  
  if (script.empty()){
    // if no script provided, read from STDIN
    shell.ExecuteStream(std::cin, "standard in");
  }else {
    logstream(LOG_EMPH) << "Flying script [" << script << "]" << std::endl;
    shell.ExecuteFile(script.c_str());
  }

  return EXIT_SUCCESS;

}

int main(int argc, char **argv){

  global_logger().set_log_level(LOG_DEBUG);
  
  graphlab::command_line_options clopts("GraphLab Javascript Shell");
  std::string script;
  clopts.attach_option("script", &script, script,
                       "The javascript file.  If none is provided, then I will read "
                       "the script from STDIN");
 
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  
  // graphlab - build communication layers
  mpi_tools::init(argc, argv);
  pilot::set_clopts(clopts);

  int rc = EXIT_FAILURE;
  try {
    logstream(LOG_INFO) << "Initializing pilot ..." << std::endl;
    rc = v8_main(script);
    logstream(LOG_INFO) << "Launch completed." << std::endl;
  } catch (std::exception const & e){
    logstream(LOG_ERROR) << "Exception : " << e.what() << std::endl;
    logstream(LOG_ERROR) << "Launch failed." << std::endl;
  }
  
  // graphlab - teardown communication layers
  mpi_tools::finalize();
  logstream(LOG_INFO) << "Bye." << std::endl;
  
  return rc;

}
