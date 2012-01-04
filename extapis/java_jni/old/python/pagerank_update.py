import math

# Python implementation of pagerank

damping = 0.85

def update(scope, scheduler):
    pvertex = scope.getVertex().value
    sumval = pvertex.rank * pvertex.selfedge + sum([e.value * scope.getNeighbor(e.from).value.rank for e in scope.getInboundEdges()])
    newval = (1-damping)/scope.getNumOfVertices() + damping*sumval
    if (abs(newval-pvertex.rank)>0.00001):
        scheduler.addTaskToOutbound(scope)
    pvertex.rank = newval





update(scope, scheduler)