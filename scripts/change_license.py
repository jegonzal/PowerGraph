#!/usr/bin/python



import os
import sys

def update_source(filename, oldcopyright, copyright):
    fdata = file(filename,"r+").read()
    # If there was a previous copyright remove it
    if (oldcopyright != None):
        if (fdata.startswith(oldcopyright)):
            fdata = fdata[len(oldcopyright):]
    # If the file does not start with the new copyright 
    if not (fdata.startswith(copyright)):
        print "  updating: " + filename
        fdata = copyright + fdata
        file(filename,"w").write(fdata)

def recursive_traversal(dir,  oldcopyright, copyright):
    fns = os.listdir(dir)
    print "Processing directory: "+dir
    for fn in fns:
        fullfn = os.path.join(dir,fn)
        if (os.path.isdir(fullfn)):
            recursive_traversal(fullfn, oldcopyright, copyright)
        else:
            if (fullfn.endswith(".cpp") or fullfn.endswith(".hpp") or
                fullfn.endswith(".cxx") ):
                update_source(fullfn, oldcopyright, copyright)


oldcright = file(sys.argv[1],"r+").read()
cright = file(sys.argv[2],"r+").read()
recursive_traversal(sys.argv[3], oldcright, cright)

exit()
