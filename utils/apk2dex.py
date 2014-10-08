import os
import sys
import getopt
from zipfile import ZipFile

"""
This script needs the latest version of zipfile to work, please get one from:
https://hg.python.org/cpython/file/2.7/Lib/zipfile.py
And replace the default zipfile.py at /usr/lib/python2.7/zipfile.py the downloaded one if default is not the latest.
"""

def process_dir(dirname, outputdir, sha1=False):
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    for filename in os.listdir(dirname):
        try:
	    targetname = filename.replace(".zip", ".apk").replace(".apk", ".dex")
	    targetname = os.path.join(outputdir, targetname)
            zfile = ZipFile(os.path.join(dirname, filename))
            for name in zfile.namelist():
                (dname, filename) = os.path.split(name)
                if filename == "":
		    # directory
		    pass
	        else:
		    # file
		    if name.endswith(".dex"):
		        with open(targetname, "w") as f:
			    f.write(zfile.read(name))
	    zfile.close()
        except:
            print "Unzip file %s error" % filename 
            
def usage():
    print """
    Usage:
    -----
        python %s -i <inputdir> -o [outputdir]
    
    Valid options are:

	-h 	You are looking at this
     	-i 	Input directory
	-o 	Output directory
    """ % (sys.argv[0])

def main():    
    if len(sys.argv) < 3:
        usage()
   	sys.exit(0)
    else:
	try:
	    opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["help", "inputdir", "outputdir"])
        except getopt.GetoptError, err:
            print str(err)
	    usage()
	    sys.exit(1)

    outputdir = "dexs"
    inputdir = ""

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("-i", "--inputdir"):
            inputdir = a
        elif o in ("-o", "--outputdir"):
            outputdir = a

    if inputdir != "":
        process_dir(inputdir, outputdir)

if __name__ == "__main__":
    main()
