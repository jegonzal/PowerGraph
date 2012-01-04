

class pagerank_vertex:
    def __init__(self, value, selfedge):
        self.rank = value
        self.selfedge = selfedge

f = open(filename, "r")
lines = f.readlines()

# First line is header
header = lines[0]
lines = lines[1:]

nvertices = int(header.split(",")[1])
for i in range(0,nvertices):
   # format: first value is value, second is self edge weight
   graph.addVertex(pagerank_vertex(1.0/nvertices, 0.0))

for l in lines:
    t = l.split(",")
    i = int(t[0])-1
    j = int(t[1])-1
    # Each line ends (annoyingly) to \n
    w =  float(t[2][:-1])
    if (i != j):
        graph.addEdge(j, i, w)
    else:
        graph.getVertex(i).value.selfedge = w