from multiprocessing import Pool
import subprocess
import os

PATTERN=('.png','.jpg','.tif')
INDEX_COUNT='4'

smoothWeight='500'
timeEnd='32'
mismatchError='0.05'
imagePath='/Users/kthierbach/Documents/current/emb/refdata_smaller/bin'
flowPath='/Users/kthierbach/Documents/current/emb/refdata_smaller/bin/flows'
fromFrame=1
toFrame=209
parallelProcesses=2



def filenames(path,pattern):
	files=filter(lambda file: any(file.lower().endswith(x) for x in pattern),os.listdir(path))
	files=map(lambda file: os.path.join(path,file),files)
	return files

def register(frame):
	subprocess.call(['fluidReg','-t ' + timeEnd,'-m ' + mismatchError,'-s ' + smoothWeight,'-f' + os.path.join(flowPath,('flow%0' + INDEX_COUNT + 'i.dat') % frame), files[frame], files[frame-1]])

def execute():
	pool=Pool(parallelProcesses)
	pool.map(register,range(fromFrame,toFrame-1))

files=filenames(imagePath,PATTERN);
execute()