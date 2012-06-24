#include "cvv8.hpp"
#include "pilot.hpp"

using namespace graphlab;
using namespace v8;
namespace cv = cvv8;

// TODO: fix this -- clopts really shouldn't be static:
graphlab_options pilot::opts;

// TODO: how many distributed controls can a pilot have in his lifetime? 
pilot::pilot() : dc(), graph(dc, opts) {}
    
void pilot::ping(){ std::cout << "pong." << std::endl; }
 
// TODO: how many graphs can a pilot have in his lifetime? 
void pilot::load_graph(const std::string &path, const std::string &format){
  graph.load_format(path, format);
  graph.finalize();
}

/**
 * Takes a javascript constructor to create vertex programs.
 */
void pilot::fly(const v8::Handle<v8::Function> &function){
  using namespace v8;
  Local<Object> obj = Object::Cast(*function->Call(function, 0, NULL));
  Local<Function> jf = Function::Cast(*(obj->Get(JSTR("apply"))));
  jf->Call(obj, 0, NULL);
}

////////////////////////////// STATIC //////////////////////////////

/**
 * Adds a JS binding of the class to the given object. Throws
 * a native exception on error.
 */
void pilot::setup_bindings(v8::Handle<v8::Object> const &dest){
  cv::ClassCreator<pilot>::Instance().SetupBindings(dest);
}

/**
 * Saves command line options for this session.
 */
void pilot::set_clopts(const graphlab_options &clopts){ opts = clopts; }

namespace cvv8 {

  CVV8_TypeName_IMPL((pilot), "pilot");

  pilot *ClassCreator_Factory<pilot>::
  Create(Persistent<Object> & jsSelf, Arguments const & argv){
    typedef CtorArityDispatcher<pilotCtors> Proxy;
    pilot * b = Proxy::Call( argv );
    if(b) BMap::Insert( jsSelf, b );
    return b;
  }

  void ClassCreator_Factory<pilot>::
  Delete(pilot *obj){
    BMap::Remove(obj);
    delete obj;
  }
  
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
