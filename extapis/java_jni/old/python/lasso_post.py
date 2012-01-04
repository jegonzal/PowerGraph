
leastsqr_err = 0.0
penalty = 0.0

lamb = 0.5

for v in graph.getVertices():
    lassov = v.value
    if lassov.vtype == 0:
        penalty += lamb * abs(lassov.value)
    else:
        leastsqr_err += pow(lassov.observed - lassov.curval,2)
        

print("Objective:", penalty + leastsqr_err)