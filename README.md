malcluster
==========

Apply fuzzy hashing techniques in Android malware analysis

1) Input
   Create a folder named "samples" in current directory, copy all the target sample files you want to analysis, the file format may depend on specific fuzzy hashing algorithms.

2) Calling fuzzy hashing algorhtms and Output
   clustering.py is the interface for calling fuzzy hashing algorithms. Conventionally, a fuzzy hashing algorithm needs to define generateHash() and compareHash() methods in order to be used in for clustering analysis. 

3) Results Evaluation
   Currently, there are three methods for evalutating fuzzy hashing analysis results: precision and recall balance point for hierarchical clustering analysis; CDF graph for inter-family and intra-family distance computation; ROC curve and AUC value for different distance threshold.
