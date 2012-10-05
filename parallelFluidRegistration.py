from multiprocessing import Pool
import subprocess
import os
import sys

PATTERN=('.png','.jpg','.tif')
INDEX_COUNT='4'

method='NavierLame'
boundary="Periodic"
timeEnd='32'
mismatchError='0.00005'

#binary
pathToFluidReg='./fluidReg'

#options for NavierLame
lamemu='1'
Lambda='.25'

#options for OverDampedCurvature and OverDampedDiffusion
smoothWeight='500'
vortexWeight='0'

imagePath='/Users/kthierbach/Documents/current/emb/refdata_smaller/gaussian/images'
flowPath='/Users/kthierbach/Documents/current/emb/refdata_smaller/gaussian/flows'
fromFrame=1
toFrame=5
parallelProcesses=2



def filenames(path,pattern):
	files=filter(lambda file: any(file.lower().endswith(x) for x in pattern),os.listdir(path))
	files=map(lambda file: os.path.join(path,file),files)
	print "Found " + str(len(files)) + " images."
	if(len(files) == 0):
		sys.exit(1)
	return files

def register(frame):
	subprocess.call([pathToFluidReg, '-t ' + timeEnd,'-e ' + mismatchError,'-s ' + smoothWeight, '-w ' + vortexWeight, '-l ' + Lambda, '-m ' + lamemu, '--flow-file', os.path.join(flowPath,('flow%0' + INDEX_COUNT + 'i.dat') % (frame+1)), files[frame], files[frame-1]])

def execute():
	pool=Pool(parallelProcesses)
	pool.map(register,range(fromFrame-1,toFrame-1))

if(os.path.isfile(pathToFluidReg) == False):
	print "Binary " + pathToFluidReg + " not found."
	sys.exit(1)
if(os.path.exists(imagePath) == False):
	print "Image path " + imagePath + " does not exist."
	sys.exit(1)
if(os.path.exists(flowPath) == False):
	print "Flow path " + flowPath + " does not exist."
	sys.exit(1)

files=filenames(imagePath,PATTERN);
execute()