import os
import sys
import getopt
import hashlib
import shutil 

def filesha1(filename):
    sha1 = hashlib.sha1()
    block_size=128*sha1.block_size
    with open(filename, 'rb') as f:
        while True:
            data = f.read(block_size)
            if not data:
                break
            sha1.update(data)
    return sha1.hexdigest()

def process_dir(dirname, outputdir, sha1=False):
    filelist = {}
    for root, dirs, files in os.walk(dirname):
        for f in files:
            targetname = root.split("/")[1] + "-" + f
            fname = os.path.join(root,f)
            filelist[fname] = os.path.join(outputdir, targetname)
    if sha1:
        for fname in filelist:
            targetname = filelist[fname].split("-")[0] + filesha1(fname)
    
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    for fname in filelist:
        try:
            shutil.copyfile(fname, filelist[fname]) 
        except KeyboardInterrupt:
            sys.exit(0)

def usage():
    print """
    Usage:
    -----
        python %s -i <inputdir> -o [outputdir] -s
    
    Valid options are:

	-h 	You are looking at this
     	-i 	Input directory
	-o 	Output directory
	-s	Include sha1 hash for each file
    """ % (sys.argv[0])

def main():    
    if len(sys.argv) < 3:
        usage()
   	sys.exit(0)
    else:
	try:
	    opts, args = getopt.getopt(sys.argv[1:], "hsi:o:", ["help", "sha1", "inputdir", "outputdir"])
        except getopt.GetoptError, err:
            print str(err)
	    usage()
	    sys.exit(1)

    outputdir = "samples"
    inputdir = ""
    computesha1 = False

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("-i", "--inputdir"):
            inputdir = a
        elif o in ("-o", "--outputdir"):
            outputdir = a
        elif o in ("-s", "--sha1"):
            computesha1 = True

    if inputdir != "":
        process_dir(inputdir, outputdir, computesha1)

if __name__ == "__main__":
    main()
