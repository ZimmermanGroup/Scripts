import numpy as np
import sys
import os.path

pobf_file = os.path.exists("xx00.pobf")

molfile = sys.argv[1]

noH = open("noH.%s" %molfile)
noH_XYZ = np.genfromtxt( noH, dtype = None, autostrip = True, skip_header = 2, usecols = (1, 2, 3))

atomsfile = open("%s" %molfile)
atoms = np.genfromtxt( atomsfile, autostrip = True, skip_header = 2, usecols = (0), dtype = str)

wholefile = open("%s" %molfile)
XYZ = np.genfromtxt( wholefile, dtype = None, autostrip = True, skip_header = 2, usecols = (1, 2, 3))

avgX = np.mean(noH_XYZ[:,0])
avgY = np.mean(noH_XYZ[:,1])
avgZ = np.mean(noH_XYZ[:,2])

centroid = np.array([avgX, avgY, avgZ])

XYZ = XYZ - centroid
noH_XYZ = noH_XYZ - centroid		# operates on all rows of XYZ

def rotation1():
	X = noH_XYZ[:,0]
	Y = noH_XYZ[:,1]
	XY = np.array([X, Y])
	XY = np.transpose(XY)
	d_matrix = []
	max_r = c1 = c2 = 0
	print pobf_file
	if pobf_file == True:
		num1 = open("xx00.pobf", "r").read().splitlines()[0]
		num1 = int(num1)
		num2 = open("xx00.pobf", "r").read().splitlines()[1]
		num2 = int(num2)
	else:
		for item1 in XY:
			c1 += 1
			c2 = 0
			for item2 in XY:
				c2 += 1
				d_matrix += [float(np.sqrt(np.sum((item1-item2)**2)))]
				if float(np.sqrt(np.sum((item1-item2)**2))) >= max_r:
					max_r = float(np.sqrt(np.sum((item1-item2)**2)))
					atom1 = item1
					atom2 = item2
					num1 = c1
					num2 = c2
				else:
					pass
	v1xyz = noH_XYZ[[num1-1]]
	print "v1xyz: ", v1xyz
	noH_trans = noH_XYZ - v1xyz
	v2xyz = noH_trans[[num2-1]]
	v2xyz = v2xyz/np.linalg.norm(v2xyz)
	print "v2xyz: ", v2xyz
	v2x, v2y, v2z = v2xyz.ravel()
	print "v2x: ", v2x
	print "v2y: ", v2y
	print "v2z: ", v2z
	if v2y == 0:
		theta = np.pi/2
	else:
		theta = np.arctan(-v2z/v2y)
	cost = np.cos(theta)
	sint = np.sin(theta)
	alpha = np.arctan(-(cost*v2y-sint*v2z)/v2x)
	cosa = np.cos(alpha)
	sina = np.sin(alpha)
	R1 = np.array([[cosa, -sina*cost, sina*sint],[sina, cosa*cost, -cosa*sint],[0, sint, cost]])
	sx, sy, sz = np.dot(R1, v2xyz.transpose()).ravel()
	if sx < 0:
		print "R1 flip"
		ntheta = np.pi 
		cosnt = np.cos(ntheta)
		sinnt = np.sin(ntheta)
		R1 = np.dot( np.array([[cosnt, -sinnt, 0],[sinnt, cosnt, 0],[0, 0, 1]]), R1)
	return R1, num1, num2, v1xyz

R1, num1, num2, v1xyz = rotation1()

XYZ = XYZ - v1xyz
noH_XYZ = noH_XYZ - v1xyz
R1_XYZ = np.dot(R1, XYZ.transpose()).transpose()
R1_noH_XYZ = np.dot(R1, noH_XYZ.transpose()).transpose()

print "num1:", num1
print "num2:", num2
print "v2xyz R1:", R1_noH_XYZ[[num2-1]]
#num1 = str(num1)
#num2 = str(num2)

name = "%s" %molfile

def tictictic():
	resnum = int(name.translate(None, 'x'))
	p_resnum = resnum - 1
	if os.path.exists( "%02d" %p_resnum + ".res") == True:
		p_resnum_file = open( "%02d" %p_resnum + ".res")
		p_resnum_matrix = np.genfromtxt( p_resnum_file, dtype = None)
		oldX = p_resnum_matrix[:,0]
		oldY = p_resnum_matrix[:,1]
		oldZ = p_resnum_matrix[:,2]
		theta = i  = j = Xres = Yres = Zres = 0
		min_Res = 1000000000000000000000000000   #just a big number that will be replaced
		for t in range(720):
			theta = (t * 2 * np.pi) / 720
			sint = np.sin(theta)
			cost = np.cos(theta) 
			RX = np.array([[1,0,0],[0, cost, -sint],[0, sint, cost]])
			X = np.dot(RX, R1_noH_XYZ.transpose()).transpose()[:,0]
			Y = np.dot(RX, R1_noH_XYZ.transpose()).transpose()[:,1]
			Z = np.dot(RX, R1_noH_XYZ.transpose()).transpose()[:,2]
			Xres = Yres = Zres = i = j = k = 0
			for item1 in X:
				Xres += (X[i] - oldX[i]) ** 2
				i += 1
			for item2 in Y:
				Yres += (Y[j] - oldY[j]) ** 2
				j += 1
			for item3 in Z:
				Zres += (Z[k] - oldZ[k]) ** 2
				k += 1
			TotalRes = Xres + Yres + Zres
			if min_Res >= TotalRes:
				min_Res = TotalRes
				min_theta = theta
			#print "TotalRes:", TotalRes
			#print "t:", t
			#print "theta:", theta
		sinmt = np.sin(min_theta)
		cosmt = np.cos(min_theta)
		print "min_Res:", min_Res
		print "min_theta:", min_theta
		Rxfinal = np.array([[ 1, 0, 0],[ 0, cosmt, -sinmt],[ 0, sinmt, cosmt]])
	else:
		Rxfinal = np.array([[1,0,0],[0,1,0],[0,0,1]])
	tic = np.dot(Rxfinal, R1_noH_XYZ.transpose()).transpose()
	np.savetxt("%02d" %resnum + ".res", tic)
	return Rxfinal

R2 = tictictic()

R2_XYZ = np.dot(R2, R1_XYZ.transpose()).transpose()
R2_noH_XYZ = np.dot(R2, R1_noH_XYZ.transpose()).transpose()

R2_XYZ = R2_XYZ - R2_XYZ[0]
R2_noH_XYZ = R2_noH_XYZ - R2_noH_XYZ[0]

print "v2xyz R2:", R2_noH_XYZ[[num2-1]]

newXYZ = np.char.mod("%.7f", R2_XYZ)
finalXYZ = np.column_stack((atoms, newXYZ))

# create a new file with the atoms and rotated/translated xyz coordinates

header_file = open(name)
for i, line in enumerate(header_file):
	if i == 0:
		top = line
	if i == 1:
		bottom = line
	else:
		pass

num1 = str(num1)
num2 = str(num2)

if pobf_file == False:
	with open(molfile + ".pobf", "w") as long_axis_file:
		long_axis_file.write(num1 + "\n")
		long_axis_file.write(num2)

with open(name + ".flat", "w") as f2:
	np.savetxt(f2, finalXYZ, fmt ='%.9s')
with open(name + ".flat", "r+") as f2:	
	old = f2.read()					# make temp file of original file
	f2.seek(0)						# rewind the file
	f2.write(top  + bottom + old)  # add header and empty line

print "---------------------------------------"
#print "Plane of Best Fit:   ",fitPlaneSVD(noH_XYZ)
