#!/usr/bin/env python

"""
Check details of clustering results
"""

def main():
    clusters = {}
    with open("eval/cluster_results.txt") as f:
	currentfamily = ""
	for line in f:
	    if line.rstrip() == "":
		continue
	    if line.startswith("C"):
		currentfamily = line.rstrip().replace(":", "")
		clusters[currentfamily] = []
	    else:
		clusters[currentfamily].append(line.rstrip())	    
    for i in clusters:
	output = ""
	results = {} 
	for mal in clusters[i]:
	    if mal.split("-")[0] not in results:
		results[mal.split("-")[0]] = [mal]
	    else:
		results[mal.split("-")[0]].append(mal)
	for family in results:
	    output += family + ":" + str(len(results[family])) + "  "
        print "Cluster ", str(i), ":\t", output
    print "Done", len(clusters)

if __name__ == "__main__":
    main()
