import cPickle;

usermod = None;

def initUserModule(name):
	global usermod;
	usermod = __import__(name);
        return "edgeDataClass" in dir(usermod);

def newVertex():
	return usermod.vertexDataClass();

def loadVertex(vertexWrap):
        return cPickle.loads(vertexWrap);

def storeVertex(vertex):
	return cPickle.dumps(vertex);

def newEdge():
	return usermod.edgeDataClass();

def loadEdge(edgeWrap):
	return cPickle.loads(edgeWrap);

def storeEdge(edge):
	return cPickle.dumps(edge);

def newAgg():
	return usermod.aggregatorClass();

def loadAgg(aggWrap):
	return cPickle.loads(aggWrap);

def storeAgg(agg):
	return cPickle.dumps(agg);

def gatherAgg(agg1, agg2):
	agg1.merge(agg2);
	return agg1;

def transformVertex(vertex):
	return usermod.transformVertex(vertex);

def transformEdge(edge):
        return usermod.transformEdge(edge);

def saveVertex(vertex):
	return usermod.saveVertex(vertex);

def saveEdge(edge):
        return usermod.saveEdge(edge);

def gather(srcData, targetData, edgeData, numIn, numOut):
        return usermod.aggregatorClass(usermod.gather(srcData, targetData, edgeData, numIn, numOut));

def apply(targetData, agg, numIn, numOut):
        return usermod.apply(targetData, agg, numIn, numOut);

def scatter(srcData, targetData, edgeData, numIn, numOut):	
        return usermod.scatter(srcData, targetData, edgeData, numIn, numOut);
