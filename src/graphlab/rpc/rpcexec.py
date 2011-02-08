#!/usr/bin/python
import sys
import os
import string
import subprocess

"""
Usage: rpcexec -n n_to_start -f [hostsfile] [program] [options]
To start local only: rpcexec [program] [options]
"""

def escape(s):
  s = string.replace(s, '"', '\\"')
  s = string.replace(s, "'", "\\'")
  return s
#enddef



nmachines = 0
hostsfile = ''
prog = ''
opts = ''
printhelp = 0
i = 1
while(i < len(sys.argv)):
  if sys.argv[i] == '-h':
    printhelp = 1
    break
  elif sys.argv[i] == '-n':
    nmachines = int(sys.argv[i+1])
    i = i + 2
  elif sys.argv[i] == '-f':
    hostsfile = sys.argv[i+1]
    i = i + 2
  else:
    prog = sys.argv[i]
    if (len(sys.argv) > i+1):
      opts = string.join(sys.argv[(i+1):])
    #endif
    break
  #endif
#endwhile

if (printhelp):
  print
  print("Usage: rpcexec -n n_to_start -f [hostsfile] [program] [options]")
  print("To start local only: rpcexec [program] [options]")
  print
  exit(0)
#endif

if (nmachines == 0 and hostsfile == ''):
  cmd = 'env SPAWNNODES=localhost SPAWNID=0 %s %s' % (prog, opts)
  p = subprocess.Popen(cmd, shell=True)
  os.waitpid(p.pid, 0)
  exit(0)
#endif
print('Starting ' + str(nmachines) + ' machines')
print('Hosts file: ' + hostsfile)
print('Command Line to run: ' + prog + ' ' + opts)

# construct the command line
sshcmd = 'ssh -x -n -q '


# open the hosts file and read the machines
try:
  f = open(hostsfile, 'r')
except:
  print
  print("Unable to open hosts file")
  print
  exit(0)
#endtry

machines = [''] * nmachines
port = [22] * nmachines
for i in range(nmachines):
  try:
    machines[i] = string.strip(f.readline())
    colonsplit = string.split(machines[i], ':')
    if (len(colonsplit) == 2):
      machines[i] = string.strip(colonsplit[0])
      port[i] = int(colonsplit[1])
    #endif
  except:
    print
    print("Unable to read line " + str(i+1) + " of hosts file")
    print
    exit(0)
#endfor
f.close()

allmachines = '"' + string.join(machines, ',') + '"'
# now issue the ssh commands
procs = [None] * nmachines
cwd = os.getcwd()
allmachines = allmachines
prog = prog
opts = opts
for i in range(nmachines):
  if (machines[i] == "localhost" or machines[i].startswith("127.")):
    cmd = 'env SPAWNNODES=%s SPAWNID=%d %s %s' % (allmachines,i, prog, opts)
  elif (port[i] == 22):
    cmd = sshcmd + '%s "cd %s ; env SPAWNNODES=%s SPAWNID=%d %s %s"' % (machines[i], escape(cwd), escape(allmachines),i, escape(prog), escape(opts))
  else:
    cmd = sshcmd + '-oPort=%d %s "cd %s ; env SPAWNNODES=%s SPAWNID=%d %s %s"' % (port[i], machines[i], escape(cwd), escape(allmachines),i, escape(prog), escape(opts))
  #endif
  print cmd
  procs[i] = subprocess.Popen(cmd, shell=True)
#endfor

for i in range(nmachines):
  os.waitpid(procs[i].pid, 0)
#endfor
