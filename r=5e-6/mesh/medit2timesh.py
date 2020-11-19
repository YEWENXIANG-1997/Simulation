from math import sqrt
from sys import argv

def readMeditMesh (meshfile="corner.mesh"):
   infile = open(meshfile,"r")
   lines = infile.readlines()
   infile.close()
   nodelist = []
   numNodes = int(lines[4].split()[0])
   print(numNodes)
   for i in range(numNodes):
      nodelist.append([0,0.0,0.0,0.0])

   id = 0
   for line in lines[5:numNodes+5]:
      data = line.split()
      X = float(data[0])
      Y = float(data[1])
      Z = float(data[2])
      nodelist[id][0] = id
      nodelist[id][1] = X
      nodelist[id][2] = Y
      nodelist[id][3] = Z
      id += 1
   numEdge = int(lines[numNodes+6].split()[0])+2
   numTriangles = int(lines[numNodes+6+numEdge].split()[0])
   triList = []
   print(numTriangles)
   for i in range(numTriangles):
      triList.append([0,0,0,0])

   tid = 0
   for line in lines[numNodes+7+numEdge:numNodes+7+numEdge+numTriangles]:
      data = line.split()
      print(data)
      id1 = int(float(data[0]))-1
      id2 = int(float(data[1]))-1
      id3 = int(float(data[2]))-1
      triList[tid][0] = tid
      triList[tid][1] = id1
      triList[tid][2] = id2
      triList[tid][3] = id3
      tid += 1

   return (nodelist, triList)

def writeTriMesh (nodeList, triList, outfilename="trimesh.lsm"):
   outstring = "Triangle\n3D-Nodes %d\n" % (len(nodeList))
   for node in nodeList:
      outstring += "%d %d 0 %f %f %f\n" % (node[0],node[0],node[1],node[2],node[3])
   outstring += "\nTri3 %d\n" % (len(triList))
   for tri in triList:
      outstring += "%d 0 %d %d %d\n" % (tri[0],tri[1],tri[2],tri[3])

   outfile = open(outfilename,"w")
   outfile.write(outstring)
   outfile.close()
   return

def normalVector(p1,p2,p3):
   r1 = [p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]]
   r2 = [p3[0]-p1[0],p3[1]-p1[1],p3[2]-p1[2]]

   cp = [r2[1]*r1[2] - r2[2]*r1[1], r2[2]*r1[0] - r2[0]*r1[2], r2[0]*r1[1] - r2[1]*r1[0]]

   norm = sqrt(cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2])

   if (norm > 0.0):
      nV = [cp[0]/norm,cp[1]/norm,cp[2]/norm]
   else:
      nV = [0,0,0]

   return nV

def reorderNodes(nodeList, triList, inPoint = [0.5,0.5,1.0]):
   for tid in range(len(triList)):
      print(tid)    
      id1 = triList[tid][1]
      id2 = triList[tid][2]
      id3 = triList[tid][3]
      p1 = [nodeList[id1][1], nodeList[id1][2], nodeList[id1][3]]
      p2 = [nodeList[id2][1], nodeList[id2][2], nodeList[id2][3]]
      p3 = [nodeList[id3][1], nodeList[id3][2], nodeList[id3][3]]

      nV = normalVector(p1,p2,p3)

      dp = [inPoint[0]-p1[0],inPoint[1]-p1[1],inPoint[2]-p1[2]]

      dp_nV = dp[0]*nV[0] + dp[1]*nV[1] + dp[2]*nV[2]

      if (dp_nV < 0):
         triList[tid][1] = id3 
         triList[tid][3] = id1 

   return

def Medit2TriMesh (meshfile="corner.mesh", outfile="corner.lsm",inPoint = [0.5,0.5,1.0]):
   (nodeList, triList) = readMeditMesh(meshfile=meshfile)

   reorderNodes (nodeList, triList, inPoint = inPoint)

   writeTriMesh (nodeList,triList,outfilename=outfile)

   return

if __name__ == "__main__":
   if (len(argv)>=3):
      Medit2TriMesh(meshfile=argv[1],outfile=argv[2],inPoint=[0,10,0])

   else:
      Medit2TriMesh(meshfile="cylinder",outfile="cylinder.lsm",inPoint=[0,0,0.001])

# save as mesh as INRIA medit
