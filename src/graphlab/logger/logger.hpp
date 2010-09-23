/**
 * @file logger.hpp
 * Usage:
 * First include logger.hpp. To logger, use the logger() function
 * There are 2 output levels. A "soft" output level which is
 * set by calling global_logger.set_log_level(), as well as a "hard" output
 * level OUTPUTLEVEL which is set in the source code (logger.h).
 *
 * when you call "logger()" with a loglevel and if the loglevel is greater than
 * both of the output levels, the string will be written.
 * written to a logger file. Otherwise, logger() has no effect.
 *
 * The difference between the hard level and the soft level is that the
 * soft level can be changed at runtime, while the hard level optimizes away
 * logging calls at compile time.
 *
 * @author Yucheng Low (ylow)
 */
#ifndef GRAPHLAB_LOG_LOG_HPP
#define GRAPHLAB_LOG_LOG_HPP
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
/**
 * \def LOG_FATAL
    Used for fatal and probably irrecoverable conditions
 * \def LOG_ERROR
    Used for errors which are recoverable within the scope of the function
 * \def LOG_WARNING
    Logs interesting conditions which are probably not fatal
 * \def LOG_INFO
    Used for providing general useful information
 *
 */
#define LOG_FATAL 3
#define LOG_ERROR 2
#define LOG_WARNING 1
#define LOG_INFO 0

/**
 * \def OUTPUTLEVEL
 *  The minimum level to logger at
 * \def LOG_NONE
 *  OUTPUTLEVEL to LOG_NONE to disable logging
 */
#define LOG_NONE 4
#define OUTPUTLEVEL LOG_INFO

/// If set, logs to screen will be printed in color
#define COLOROUTPUT


/**
 * \def logger(lvl,fmt,...)
 *    extracts the filename, line number
 *     and function name and calls _log. It will be optimized
 *     away if LOG_NONE is set
 *     This relies on a few compiler macros. As far as I know, these
 *     macros are pretty standard among most other C++ compilers. 
 */
#if OUTPUTLEVEL == LOG_NONE
// totally disable logging
#define logger(lvl,fmt,...)
#define logbuf(lvl,fmt,...)
#define logstream
#else
#define logger(lvl,fmt,...)                 \
    (global_logger()._log(lvl,__FILE__,     \
                        __func__ ,__LINE__,fmt,##__VA_ARGS__))

#define logbuf(lvl,buf,len)                 \
    (global_logger()._log(lvl,__FILE__,     \
                        __func__ ,__LINE__,buf,len))
#define logstream(lvl)                      \
    global_logger().start_stream(lvl,__FILE__, __func__ ,__LINE__) 
#endif


/**
  logging class.
  This writes to a file, and/or the system console.
*/
class file_logger{
 public:
  /** Default constructor. By default, log_to_console is off,
      there is no logger file, and logger level is set to LOG_WARNING
  */
  file_logger();  

  ~file_logger();   /// destructor. flushes and closes the current logger file

  /** Closes the current logger file if one exists.
      if 'file' is not an empty string, it will be opened and 
      all subsequent logger output will be written into 'file'.
      Any existing content of 'file' will be cleared.
      Return true on success and false on failure.
  */
  bool set_log_file(std::string file);

  /// If consolelog is true, subsequent logger output will be written to stderr
  void set_log_to_console(bool consolelog) {
    log_to_console = consolelog;
  }

  /// Returns the current logger file.
  std::string get_log_file(void) {
    return log_file;
  }

  /// Returns true if output is being written to stderr
  bool get_log_to_console() {
    return log_to_console;
  }

  /// Returns the current logger level
  int get_log_level() {
    return log_level;
  }

  file_logger& start_stream(int lineloglevel,const char* file,const char* function, int line);
  
  template <typename T>
  file_logger& operator<<(T a) {
    if (streamactive) streambuffer << a;
    return *this;
  }

  file_logger& operator<<(const char* a) {
    if (streamactive) {
      streambuffer << a;
      if (a[strlen(a)-1] == '\n') {
        stream_flush();
      }
    }
    return *this;
  }

  file_logger& operator<<(std::ostream& (*f)(std::ostream&)){
    typedef std::ostream& (*endltype)(std::ostream&);
    if (streamactive) {
      if (endltype(f) == endltype(std::endl)) {
        streambuffer << "\n";
        stream_flush();
      }
    }
    return *this;
  }



  /** Sets the current logger level. All logging commands below the current
      logger level will not be written. */
  void set_log_level(int new_log_level) {
    log_level = new_log_level;
  }

  /**
  * logs the message if loglevel>=OUTPUTLEVEL
  * This function should not be used directly. Use logger()
  *
  * @param loglevel Type of message \see LOG_INFO LOG_WARNING LOG_ERROR LOG_FATAL
  * @param file File where the logger call originated
  * @param function Function where the logger call originated
  * @param line Line number where the logger call originated
  * @param fmt printf format string
  * @param ... The parameters that match the format string
  */
  void _log(int loglevel,const char* file,const char* function,
                int line,const char* fmt, ... );

  void _logbuf(int loglevel,const char* file,const char* function,
                int line,  const char* buf, int len);
                
  void _lograw(int loglevel, const char* buf, int len);

  void stream_flush() {
    streambuffer.flush();
    _lograw(streamloglevel,
            streambuffer.str().c_str(),
            streambuffer.str().length());
    streambuffer.str("");
  }
 private:
  std::ofstream fout;
  std::string log_file;
  std::stringstream streambuffer;
  int streamloglevel;
  bool streamactive;
  
  bool log_to_console;
  int log_level;

};


file_logger& global_logger();

#include <graphlab/logger/assertions.hpp>

#endif
