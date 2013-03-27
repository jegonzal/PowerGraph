class vertexDataClass:
	pr = 1.0;
	prDelta = 0.0;
	def __init__(self, pr_new=1.0, delta_new=0.0):
		self.pr = pr_new;
		self.prDelta = delta_new;		
			
class aggregatorClass:
	sum = 0.0;
	def __init__(self, sum_new=0.0):
		self.sum = sum_new;
	def merge(self, x):
		self.sum += x.sum;

def parseEdge(file, line):		
	s = line.split(" ");
	srcId = int(s[0]);
	destId = int(s[1]);
	return (srcId, destId, None);
	
def transformVertex(vertex):
	return vertexDataClass();

def saveVertex(vertex):
	return str(vertex.pr)+":"+str(vertex.prDelta);

def gather(srcData, targetData, edgeData, numIn, numOut):
	return 0.85*srcData.pr/numOut;

def apply(targetData, aggInst, numIn, numOut):
	newval = aggInst.sum+0.15;
	delta = newval-targetData.pr;
	return vertexDataClass(newval, delta);

def scatter(srcData, targetData, edgeData, numIn, numOut):
    if abs(srcData.prDelta) > 0.01:
        return (1.0, None, None);
    else:
        return (-1.0, None, None);

