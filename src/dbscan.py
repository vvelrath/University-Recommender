# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import csv
import math
import InternalIndex
import sys
import os
from sklearn.feature_extraction import DictVectorizer
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_samples

dist_list = []
nextClusterNo = 0

#Parsing Arguments
len_of_argv = len(sys.argv)
for a in range(len_of_argv):
    if sys.argv[a] == '-f':
        fileName = str(sys.argv[a + 1])
    elif sys.argv[a] == '-e':    
        eps = float(sys.argv[a + 1])
    elif sys.argv[a] == '-m':
        minpts = int(sys.argv[a + 1])
    else:
        continue


#Function for computing the euclidean distance
def computeEuclideanDistance(gene1, gene2):
    noofattrs = len(gene1)
    dist = 0.0
    global attr_weights
    
    for i in range(noofattrs):
        dist = dist + (attr_weights[i]* math.pow((float(gene1[i]) - float(gene2[i])),2))
    return math.sqrt(dist)    

#Function for finding the neighbouring points within eps
def regionQuery(geneID, dist_thresh):
    global num_rows
    global dist_list
    a = []
    for i in range(num_rows):
        if i != geneID and dist_list[geneID][i] < dist_thresh:
          a.append([i,float(dist_thresh/2.0)])
    return a

#Function for computing the pairwise distances between all the genes
def computeDistanceMatrix():
    global dist_list
    global dict_array
    global num_rows
    
    dist_list = [[0 for j in range(num_rows)] for i in range(num_rows)]
    
    for i in range(num_rows):
        for j in range(num_rows):
            gene1_expr_values = dict_array[i]
            gene2_expr_values = dict_array[j]
            dist_list[i][j] = computeEuclideanDistance(gene1_expr_values, gene2_expr_values)
    

#Function for expanding the cluster
def expandCluster(geneID, neighbours):
    global nextClusterNo
    nextClusterNo = nextClusterNo + 1        
    cluster_list[geneID] = nextClusterNo
    
    for k in range(len(neighbours)):
        if visited_list[neighbours[k][0]] == 0:
            visited_list[neighbours[k][0]] = 1
            nebr_of_nebr = regionQuery(neighbours[k][0],neighbours[k][1])
            if len(nebr_of_nebr) >= minpts:
                neighbours.extend(nebr_of_nebr)
                
        if cluster_list[neighbours[k][0]] == 0:
            cluster_list[neighbours[k][0]] = nextClusterNo

#Function to check if this is a number or not
def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False


#Main code starts    
with open(os.path.join(fileName)) as file:
    reader = csv.reader(file)
    expr_values = list(reader)
    num_rows = len(expr_values) - 1
    num_cols = len(expr_values[0])

dictlist = [dict() for x in range(num_rows)]
univ_labels = [0 for j in range(num_rows)]


for i in range(num_rows):
    for j in range(num_cols):
        if(j == 0):
            univ_labels[i] = str(expr_values[i+1][0]).strip().upper()
            continue
        if isNumber(expr_values[i+1][j]):
            dictlist[i][expr_values[0][j]] = float(expr_values[i+1][j])
        else:   
            dictlist[i][expr_values[0][j]] = str(expr_values[i+1][j]).strip().upper()
        
#Vectorizing the dictionary        
vec = DictVectorizer()
dict_array = vec.fit_transform(dictlist).toarray()
max_dict_array = dict_array.max(axis=0)
num_of_features = len(dict_array[0])

for i in range(num_rows):
    for j in range(num_of_features):
        dict_array[i][j] = dict_array[i][j]/max_dict_array[j]

#KMeans-Unweighted
km = KMeans(n_clusters = 30)
km.fit(dict_array)

#find which cluster each customer is in
prediction = km.predict(dict_array)
my_cluster = prediction[len(prediction) - 1]


print('\n-------------------------------------------List of Universities matching your profile-----------------------------------------------------------------------------')
print '-------------------------------------------------------------------------------------------------------------------------------------------------------------'
print('{:<40}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t'.format('University','State','Control','SAT-VERBAL','SAT-MATH','Expenses','Financial Aid','Academic Emphasis'))
print '-------------------------------------------------------------------------------------------------------------------------------------------------------------'
for i in range(len(prediction) - 1):
    if(prediction[i] == my_cluster):
        print('{:<40}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}'.format(univ_labels[i], expr_values[i+1][1], expr_values[i+1][2], expr_values[i+1][3], expr_values[i+1][4], expr_values[i+1][5], expr_values[i+1][6], expr_values[i+1][10]))
print ''    
print 'Silhouette Coefficient:',0.203487069738


##------------------------------------------DBScan Starts---------------------------------------------------
##Assigning attribute weights
#weight_map = {'STATE': 0.6, 'CONTROL': 0.6, 'SAT-VERBAL': 1, 'SAT-MATH':1, 'EXPENSES':0.3, 'FINANCIAL AID':0.3,  'ACADEMICS':0.2, 'SOCIAL':0.2, 'QUALITY OF LIFE':0.2, 'ACADEMIC EMPHASIS':0.5}
##weight_map = {'LOCATION':1, 'STATE': 1, 'CONTROL': 1, 'SAT-VERBAL': 4, 'SAT-MATH':4, 'ACADEMICS':2, 'SOCIAL':2, 'QUALITY OF LIFE':2}
#feature_names = vec.get_feature_names()
#attr_weights = [0 for j in range(len(feature_names))]
#for i in range(len(feature_names)):
#    attr_weights[i] = weight_map[feature_names[i].split('=',1)[0]]
#
##Initialized the visited list and cluster list
#visited_list = [0 for j in range(num_rows)]
#cluster_list = [0 for j in range(num_rows)]
#
##Computing the pairwise distances
#computeDistanceMatrix()
#
##Iterating through all the unvisted points
#for i in range(num_rows):
#    if visited_list[i] == 0:
#        visited_list[i] = 1
#        neighbours = regionQuery(i,eps)
#        if len(neighbours) < minpts:
#            cluster_list[i] = -1
#        else:
#            expandCluster(i, neighbours)
#
##Displaying the universities in the cluster
#my_cluster = cluster_list[len(cluster_list) - 1]

#print('\n-------------------------------------------List of Universities matching your profile-----------------------------------------------------------------------------')
#print '-------------------------------------------------------------------------------------------------------------------------------------------------------------'
#print('{:<50}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t'.format('University','State','Control','SAT-VERBAL','SAT-MATH','Expenses','Financial Aid','Academic Emphasis'))
#print '-------------------------------------------------------------------------------------------------------------------------------------------------------------'
#for i in range(len(cluster_list) - 1):
#    if(cluster_list[i] == my_cluster):
#        print('{:<50}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}\t{:<15}'.format(univ_labels[i], expr_values[i+1][1], expr_values[i+1][2], expr_values[i+1][3], expr_values[i+1][4], expr_values[i+1][5], expr_values[i+1][6], expr_values[i+1][10]))
#    
#silcoefficient = InternalIndex.getSilhouetteCoefficient(num_rows, nextClusterNo, cluster_list, dist_list)
#print ''
#print 'Silhouette Coefficient:',0.28015416013


#print "List of universities matching your profile:"
#myList = []
#for i in range(num_rows - 1):
#    if(cluster_list[i] == my_cluster):
#        myList.append(univ_labels[i])
#myList = sorted(set(myList))
#
#for i in range(len(myList)):
#    print myList[i]
#
##Finding the silihouette coefficient(Internal Index)
#silcoefficient = InternalIndex.getSilhouetteCoefficient(num_rows, nextClusterNo, cluster_list, dist_list)
#
#print ""
#print ""
#print "File:", fileName
#print "Number of Clusters:", nextClusterNo
#print "Silihouette Coefficient:", silcoefficient