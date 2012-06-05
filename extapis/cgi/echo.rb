#!/usr/bin/env ruby

def colorize(text, color_code)
  "\e[0;#{color_code}m#{text}\e[0m"
end

def red(text)
  colorize(text, 31)
end

def green(text)
  colorize(text, 32)
end

def cyan(text)
  colorize(text, 36)
end


begin
 
  # open new terminal 
  fd = IO.sysopen("/dev/tty", "r")
  term = IO.new(fd,"r")

  thread = Thread.new(term) { |term|
    while (line = term.readpartial 4096) do
      STDOUT.write line
      STDOUT.flush
    end
  }

  while (line = STDIN.readpartial 4096) do
    break if line.match /exit/
    STDERR.write cyan(line)
  end

  thread.exit

rescue EOFError
  STDERR.write "Pipe broken"
end
