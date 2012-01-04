import math

# Python implementation of pagerank

def update(scope, scheduler):
    vertex = scope.getVertex()
    oldval = vertex.value
    newval = vertex.value * vertex.selfEdgeWeight
    newval = newval + sum([e.weight*e.value for e in scope.getInboundEdges()])
    vertex.setValue(newval)
    if (abs(newval-oldval)>0.00001):
        scheduler.addTaskToOutbound(scope)





update(scope, scheduler)