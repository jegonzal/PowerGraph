#!/usr/bin/python
import sys
import string
import subprocess

def left_child(index) : return 2 * index + 1


def right_child(index) : return 2 * index + 2


def escape_str(text) : 
    return text.replace('\\', '\\\\').replace('"', '\\"')
    

def do_send(src, dest, path ) :
    return ' ( rsync -avz ' + \
        path + '/' + src + ' ' + \
        dest + ':' + path + ' ) '



def make_str(names, a, b, index) :
    assert(left_child(index) < len(names));
    text = '( ssh ' + names[index]
    text = text + ' " hostname; ';
    text = text + \
        do_send(a, names[left_child(index)], b)
    if right_child(index) < len(names) :
         text = text + ' & ' + \
             do_send(a, names[right_child(index)], b) + \
             ' ; wait ; '
    if left_child(left_child(index)) < len(names) :
        text = text + \
            escape_str(make_str(names, a, b, left_child(index)))
    if left_child(right_child(index)) < len(names) :
        text = text + ' & ' + \
            escape_str(make_str(names, a, b, right_child(index))) + \
            ' ; wait ; '
    text = text + ' " ) '
    return text;
    
        




host_filename = sys.argv[1]


source_path = sys.argv[2]
dest_path = sys.argv[3]
textfile = open(host_filename, "r")
lines = textfile.readlines()
machines = [ x.replace('\n', '') for x in lines]

print machines 

command = 'rsync -avz ' + source_path + ' ' + \
    machines[0] + ':' + dest_path;
print command;
subprocess.call(command, shell=True)

command = make_str(machines, source_path, dest_path, 0);
print command
subprocess.call(command, shell=True)
