import numpy;

NUMLATENT = 2;
TOLERANCE = 1e-3;
LAMBDA = 0.01;
MAX_UPDATES = 10;
MAXVAL = 1e+100;
MINVAL = -1e+100;
REGNORMAL = 1;

gatherEdges = 3;  # gather and scatter on all edges
scatterEdges = 3;

indices = numpy.triu_indices(NUMLATENT);

class vertexDataClass:
	first_update = 1;
	do_gather = 1;
	num_updates = 0;
	residual = 1.0;
	factor = numpy.zeros(NUMLATENT);
	def __init__(self, num_updates = 0, residual = 1.0):
		self.first_update = 1;
		self.do_gather = 1;
		self.num_updates = num_updates;
		self.residual = residual;
		self.factor = numpy.random.rand(NUMLATENT);

class edgeDataClass:
	type = 0;  # 0 - train, 1 - validate, 2 - predict
	obs = 0.0;
	def __init__(self, type = 0, obs = 0.0):
		self.type = type;
		self.obs = obs;

class aggregatorClass:
	XtX = numpy.zeros((NUMLATENT, NUMLATENT));
	Xy = numpy.zeros(NUMLATENT);
	is_empty = 1;
	def __init__(self, X = None, y = None):
		if (X != None) & (y != None):
                        self.XtX = numpy.zeros((NUMLATENT, NUMLATENT));
                        self.Xy = numpy.zeros(NUMLATENT);
			self.XtX[indices] = (X*X[:, numpy.newaxis])[indices];
			self.Xy = X*y;
			self.is_empty = 0;
	def merge(self, other):
		if other.is_empty == 0:
			if self.is_empty == 1:
				self.XtX = other.XtX;
				self.Xy = other.Xy;
				self.is_empty = 0;
			else:
				self.XtX[indices] += other.XtX[indices];
				self.Xy += other.Xy;

def parseEdge(file, line):
	type = 0;
	if ".train" in file:
		type = 0;
	elif ".validate" in file:
		type = 1;
	elif ".predict" in file:
		type = 2;
		
	s = line.split();	
	srcId = int(s[0]);
	destId = int(s[1]);
	edgeVal = edgeDataClass(type, float(s[2]));
	return (2*srcId, 2*destId+1, edgeVal);	

def transformVertex(vertex):
    return vertexDataClass();

def saveVertex(vertexData):
    return str(vertexData.factor);

def saveEdge(srcData, targetData, edgeData):
    if edgeData.type == 2:
        prediction = numpy.dot(srcData.factor, targetData.factor);
        return str(prediction);

def gather(srcData, targetData, edgeData, numIn, numOut):
    if edgeData.type == 0:
        return aggregatorClass(srcData.factor, edgeData.obs);
    else:
        return aggregatorClass();

def apply(targetData, aggInst, numIn, numOut):
    if (targetData.first_update == 1) & (numOut == 0):
	targetData.first_update = 0;
	targetData.do_gather = 0;
        return targetData;
    else:
   	targetData.do_gather = 1;
    
    if aggInst.is_empty == 1:
        targetData.residual = 0.0;
        targetData.num_updates = 0;
        return targetData;

    regularization = LAMBDA;
    if REGNORMAL == 1:
        regularization *= numOut;
    for i in range(0, aggInst.XtX.shape[0]):
        aggInst.XtX[i][i] += regularization;

    old_factor = targetData.factor;
    targetData.factor = numpy.linalg.solve(aggInst.XtX, aggInst.Xy);
    targetData.residual = numpy.sum(numpy.absolute(targetData.factor-old_factor)) / aggInst.XtX.shape[0];
    targetData.num_updates += 1;

    return targetData;

def scatter(srcData, targetData, edgeData, numIn, numOut):
    if srcData.do_gather == 0:
        return (-1.0, None, None);
    
    if edgeData.type == 0:
        pred = numpy.dot(srcData.factor, targetData.factor);
        error = abs(edgeData.obs-pred);
        priority = error*srcData.residual;       
        if (priority > TOLERANCE) & (targetData.num_updates < MAX_UPDATES):
            return (priority, None, None);
        else:
            return (-1.0, None, None);

