
extern "C" {
  #include <hdfs.h>
}

int main(int argc, char **argv) {

  hdfsFS fs = hdfsConnect("default", 0);
  const char* writePath = "/tmp/testfile.txt";
  hdfsFile writeFile = hdfsOpenFile(fs, writePath, O_WRONLY|O_CREAT, 0, 0, 0);
  if(!writeFile) {
    fprintf(stderr, "Failed to open %s for writing!\n", writePath);
    exit(-1);
  }
  char* buffer = "Hello, World!";
  tSize num_written_bytes = hdfsWrite(fs, writeFile, (void*)buffer, strlen(buffer)+1);
  if (hdfsFlush(fs, writeFile)) {
    fprintf(stderr, "Failed to 'flush' %s\n", writePath); 
    exit(-1);
  }
  hdfsCloseFile(fs, writeFile);
}
