#include <logger/logger.hpp>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>

file_logger& global_logger() {
  static file_logger l;
  return l;
}

const char* messages[] = {  "INFO   ",
                            "WARNING",
                            "ERROR  ",
                            "FATAL  "};


file_logger::file_logger() {
  log_file = "";
  log_to_console = true;
  log_level = LOG_WARNING;
}

file_logger::~file_logger() {
  if (fout.good()) {
    fout.flush();
    fout.close();
  }
}

bool file_logger::set_log_file(std::string file) {
  // close the file if it is open
  if (fout.good()) {
    fout.flush();
    fout.close();
    log_file = "";
  }
  // if file is not an empty string, open the new file
  if (file.length() > 0) {
    fout.open(file.c_str());
    if (fout.fail()) return false;
    log_file = file;
  }
  return true;
}



#define RESET   0
#define BRIGHT    1
#define DIM   2
#define UNDERLINE   3
#define BLINK   4
#define REVERSE   7
#define HIDDEN    8

#define BLACK     0
#define RED   1
#define GREEN   2
#define YELLOW    3
#define BLUE    4
#define MAGENTA   5
#define CYAN    6
#define WHITE   7

void textcolor(FILE* handle, int attr, int fg)
{
  char command[13];
  /* Command is the control command to the terminal */
  sprintf(command, "%c[%d;%dm", 0x1B, attr, fg + 30);
  fprintf(handle, "%s", command);
}

void reset_color(FILE* handle)
{
  char command[20];
  /* Command is the control command to the terminal */
  sprintf(command, "%c[0m", 0x1B);
  fprintf(handle, "%s", command);
}



void file_logger::_log(int lineloglevel,const char* file,const char* function,
                int line,const char* fmt, ... ){
  // if the logger level fits
  if (lineloglevel >= 0 && lineloglevel <= 3 && lineloglevel >= log_level){
    // get just the filename. this line found on a forum on line.
    // claims to be from google.
    file = ((strrchr(file, '/') ? : file- 1) + 1);

    // length of the 'head' of the string
    size_t headerlen = snprintf(NULL,0,"%s:%s(%s:%d): ",
                          messages[lineloglevel],file,function,line);

    // get the length of the actual logger
    va_list ap;
    va_start(ap,fmt);
    size_t restlen = vsnprintf(NULL,0,fmt,ap);
    // allocate the string
    size_t totallen = headerlen + restlen;
    char *str = new char[totallen + 3];

    // write the actual header
    int byteswritten = snprintf(str,totallen + 1,"%s:%s(%s:%d): ",
                          messages[lineloglevel],file,function,line);
    // write the actual logger
    va_start(ap,fmt);
    byteswritten += vsnprintf(str + byteswritten,totallen - byteswritten + 1,fmt,ap);
    snprintf(str + byteswritten,totallen - byteswritten + 3,"\n");
    // write the output
    if (fout.good()) {
      fout << str;;
    }
    if (log_to_console) {
      #ifdef COLOROUTPUT
        if (lineloglevel == LOG_FATAL) {
          textcolor(stderr, BRIGHT, RED);
        }
        else if (lineloglevel == LOG_ERROR) {
          textcolor(stderr, BRIGHT, RED);
        }
        else if (lineloglevel == LOG_WARNING) {
          textcolor(stderr, BRIGHT, GREEN);
        }
      #endif
      std::cerr << str;;
      #ifdef COLOROUTPUT
        reset_color(stderr);
      #endif
    }
    
    delete [] str;
  }
}


void file_logger::_logbuf(int lineloglevel,const char* file,const char* function,
                int line,const char* buf, int len) {
  // if the logger level fits
  if (lineloglevel >= 0 && lineloglevel <= 3 && lineloglevel >= log_level){
    // get just the filename. this line found on a forum on line.
    // claims to be from google.
    file = ((strrchr(file, '/') ? : file- 1) + 1);

    // length of the 'head' of the string
    size_t headerlen = snprintf(NULL,0,"%s:%s(%s:%d): ",
                          messages[lineloglevel],file,function,line);

    if (headerlen> 2047) {
      std::cerr << "Header length exceed buffer length!";
    }
    else {
      char str[2048];
      const char *newline="\n";
      // write the actual header
      int byteswritten = snprintf(str,2047,"%s:%s(%s:%d): ",
                            messages[lineloglevel],file,function,line);
      _lograw(lineloglevel,str, byteswritten);
      _lograw(lineloglevel,buf, len);
      _lograw(lineloglevel,newline, strlen(newline));
    }
  }
}

void file_logger::_lograw(int lineloglevel, const char* buf, int len) {
  if (fout.good()) {
    fout.write(buf,len);
  }
  if (log_to_console) {
    #ifdef COLOROUTPUT
      if (lineloglevel == LOG_FATAL) {
        textcolor(stderr, BRIGHT, RED);
      }
      else if (lineloglevel == LOG_ERROR) {
        textcolor(stderr, BRIGHT, RED);
      }
      else if (lineloglevel == LOG_WARNING) {
        textcolor(stderr, BRIGHT, GREEN);
      }
    #endif
    std::cerr.write(buf,len);
    #ifdef COLOROUTPUT
      reset_color(stderr);
    #endif
  }
}

file_logger& file_logger::start_stream(int lineloglevel,const char* file,const char* function, int line) {
  if (streambuffer.str().length() > 0) {
    stream_flush();
  }
  if (lineloglevel >= 0 && lineloglevel <= 3 && lineloglevel >= log_level){
    streambuffer << messages[lineloglevel] << ":" << file
                << "(" << function << ":" <<line<<"): ";
    streamactive = true;
    streamloglevel = lineloglevel;
  }
  else {
    streamactive = false;
  }
  return *this;
}

