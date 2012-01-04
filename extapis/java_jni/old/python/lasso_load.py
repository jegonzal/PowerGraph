#
# Solve Lasso:  min ||Ax-y||_2^2 + \lambda ||x||_1  
#
# We present Lasso as a bipartite graph. On the left side, we have
# variables x_i (predictors) and on the right side the current estimates
# for y_i = (Ax)_i. Sides are connected by edges weighted by A_ij
#


#
# Rights side of the graph. Estimate for y_i. We store
# the actual y_i as well (observed), in order to quickly
# compute prediction error.
#
class lasso_estimate_vertex:
    def __init__(self, value, observed):
        self.curval = value
        self.lastval = value
        self.observed = observed
        self.vtype = 1
        
class lasso_variable_vertex:
    def __init__(self, value):
        self.value = value
        self.covar = 0
        self.Ay = 0
        self.initialized = False
        self.vtype = 0

     

f = open(filename, "r")
lines = f.readlines()

header = lines[0].split(",")
assert(header[0] == "y")

# First read y-values
ny = int(header[1])

for i in range(1,ny+1):
    val = float(lines[i])
    graph.addVertex(lasso_estimate_vertex(0.0, val))
    
# Remove the first part
lines = lines[ny+1:]

# Read edges
header = lines[0].split(",")
print(header)
assert(header[0] == "A")

n = int(header[1])
nx = int(header[2])

# Create variables
for i in range(0,nx):
    graph.addVertex(lasso_variable_vertex(0.0))

for i in range(1,n+1):
    ln = lines[i].split(",")
    idx = int(ln[0])-1
    val = float(ln[1])
    row = idx%ny
    col = idx/ny
    # Create edge between variable x_col and y_row
    graph.addEdge(col+ny, row, val)
    assert(col+ny>row)
    #print(idx, col, row)
    
print("Data loaded")