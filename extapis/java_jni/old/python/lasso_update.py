import math

lamb = 0.5

# Python implementation of Lasso Shooting algorithm.
# min ||Ax-y||_2^2  + lambda ||x||_1


def update(scope, scheduler):
    # Of class lasso_variable_vertex or lasso_estimate_vertex
    lassov = scope.getVertex().value
    
    if (lassov.vtype == 0):
        if lassov.initialized == False:
            # Initialize covariance
            lassov.covar = 2.0*sum([e.value*e.value for e in scope.getOutboundEdges()])
            # Initialize (Ay)_i
            lassov.Ay = 2.0*sum([e.value * scope.getNeighbor(e.to).value.observed  for e in scope.getOutboundEdges()])
            lassov.initialized = True
        
        # Compute (Ax)_i
        curest = sum([e.value * scope.getNeighbor(e.to).value.curval  for e in scope.getOutboundEdges()])
        newval = soft_threshold(lamb, curest*2 - lassov.covar*lassov.value - lassov.Ay)/lassov.covar
        
       # if (newval == 0.0):
           # print("zero!")
        #if (scope.getVertex().getId() % 100 == 0):
        #    print(scope.getVertex().getId(), lassov.value, newval, curest*2 - lassov.covar*lassov.value - lassov.Ay, lassov.covar, lassov.Ay)
        
        if newval != lassov.value:
            delta = newval-lassov.value
            lassov.value = newval
            for e in scope.getOutboundEdges():
                scope.getNeighbor(e.to).value.curval += delta * e.value
    

def soft_threshold (lamb, x):
    if (x > lamb):
        return (lamb-x) 
    elif (x < lamb):
        return (-lamb-x)
    else:
        return 0.0
update(scope, scheduler)