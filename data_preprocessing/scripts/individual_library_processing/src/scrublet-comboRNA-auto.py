
#Reyka Jayasinghe modded by Austin
#reyka@wustl.edu
#Edited on 2021-08-25
#Last Edit: user input
#https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb

from sklearn import cluster
import scrublet as scr
#If error in identifying above rerun `pip install scrublet` should fix all issues
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import gzip
import pandas as pd
import random
from sklearn.cluster import KMeans 
import sys, getopt

def main(argv):
    sample = ''
    cellranger = ''
    cutoff = 0.2
    try:
        opts, args = getopt.getopt(argv,"hs:c:u:",["sname=","celloutput=", "cutoff="])
    except getopt.GetoptError:
        print ('2_Multiplet.py -s <sample> -c <cellranger> -u <cutoff>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('2_Multiplet.py -s <sample> -c <cellranger> -u <cutoff>')
            print('python 2_Multiplet.py -s RM024R1-XBn1_1 -c /diskmnt/Datasets/mmy_scratch/MetNet/Brain/cellranger/RM024R1-XBn1_1/RM024R1-XBn1_1/ -u 0.2')
            sys.exit()
        elif opt in ("-s", "--sname"):
            sample = arg
        elif opt in ("-c", "--celloutput"):
            cellranger = arg
        elif opt in ("-u", "--cutoff"):
            cutoff = arg	
    print ('Sample is "', sample)
    print ('Cellranger Directory is "', cellranger)
    print('Doublet cutoff is "', cutoff)


if __name__ == "__main__":
     main(sys.argv[1:])

#Define input arguments
cellranger=sys.argv[4]
sample=sys.argv[2]
cutoff = sys.argv[6]


#input_dir = '/diskmnt/Datasets/mmy_scratch/MetNet/Brain/cellranger/RM024R1-XBn1_1/RM024R1-XBn1_1/outs/filtered_feature_bc_matrix/'
#Only works on filtered_feature matrix - raw matrix takes too much memory
print("loading starting matrix...")
input_dir=cellranger+"/outs/filtered_feature_bc_matrix/"
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
barcodes=pd.read_csv(input_dir+'barcodes.tsv.gz',compression='gzip',header=None)
#keep gene expression features only
features=pd.read_csv(input_dir+'features.tsv.gz',compression='gzip',header=None, delimiter='\t')
index = np.array(features[features[2]=="Gene Expression"].index)
counts_matrix = counts_matrix[:,index]
#print(counts_matrix)
number_of_iterations = 10
print("initializing data sink...")
#initializing simulated doublets storage dataframe
index = range((2*len(barcodes)))
columns = []
for i in range(1,(number_of_iterations+1)):
    columns.append("test"+str(i))
sim_df = pd.DataFrame(index=index, columns=columns)
#initialize list for storing Kmeans determined cutoff
cutoff_list = []
#function for determining cutoff that should be used via kmeans clustering
def kmeans_cutoff(sim_doublets):
    cluster_df = pd.DataFrame(sim_doublets)
    km = KMeans(n_clusters=2, init='k-means++',n_init=10,max_iter=10000)
    sim_fit = km.fit_predict(cluster_df)
    cluster1_max = max(cluster_df[sim_fit==0][0])
    cluster2_max = max(cluster_df[sim_fit==1][0])
    new_cutoff = min(cluster1_max, cluster2_max)
    return new_cutoff

##########
#loop passes 1 through 9 (inclusive) which
##########
print("starting doublet prediction loop for determining cutoff")
for i in range(1,number_of_iterations):
    print("Loop "+str(i))
    seed = random.randint(0,1000000)
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.15, random_state=seed) #provide alternative random state seed for all but final loop pass.
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3,log_transform=True,
                                                            #mean_center=True,
                                                            #normalize_variance=True,
                                                            min_gene_variability_pctl=85, n_prin_comps=30)
    scrub.call_doublets(threshold=float(cutoff))
    sim_doublets = np.asarray(scrub.doublet_scores_sim_)
    sim_df['test'+str(i)]= sim_doublets
    boundary = kmeans_cutoff(sim_doublets)
    print("The boundary between clusters was "+str(boundary))
    cutoff_list.append(boundary)

##########
#10th pass through the loop done separetely to allow for recording of simulated doublets that are plotted in the final object
##########
print("Generating final doublet prediction with KMeans determined cutoff")
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.15) #initialize scrublet object with default seed of 0. This allows for consistent replication of final plots.
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3,log_transform=True,
                                                        #mean_center=True,
                                                        #normalize_variance=True,
                                                        min_gene_variability_pctl=85, n_prin_comps=30) #initial doublet prediction
scrub.call_doublets(threshold=float(cutoff)) #gives threshold
sim_doublets = np.asarray(scrub.doublet_scores_sim_)
sim_df['test'+str(number_of_iterations)]= sim_doublets
boundary = kmeans_cutoff(sim_doublets)
print("The boundary between clusters was "+str(boundary))
cutoff_list.append(boundary)
print("Generating output files and plots")
#writing simulated doublets dataframe to csv that can be manually reviewed if need be
output_simulated_doublets= sample+"_scrublet_simulated_doublet_scores_table.csv"
sim_df.to_csv(output_simulated_doublets, sep = ",", index=False)
final_cutoff = np.mean(cutoff_list)
scrub.call_doublets(threshold=float(final_cutoff)) #gives threshold
out_file = open(sample+"_scrublet_cutoff_.txt", "w")
#print(cutoff_list)
for value in cutoff_list:
    out_file.write(str(value)+"\n")
out_file.write("The final cutoff used was: "+str(final_cutoff))
out_file.close()
# scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.15) #initialize scrublet object with default seed of 0. This allows for consistent replication of final plots and inclusion of final cutoff determined from kmeans.
# doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3,log_transform=True,
#                                                         #mean_center=True,
#                                                         #normalize_variance=True,
#                                                         min_gene_variability_pctl=85, n_prin_comps=30) #initial doublet prediction
# scrub.call_doublets(threshold=float(np.mean(cutoff_list))) #gives threshold
# out_file = open(sample+"_scrublet_output_table.csv")

# centroids = np.zeros((2,))
# for i in range(0, 2):
#     print(i)
#     index=random.randint(0, len(sim_doublets))
#     print(index)
#     centroids[i,]=index

# classifications = np.zeros((sim_doublets.shape[0],), dtype=int)
# def assignPointsToCentroids():
#     for i in range(sim_doublets.shape[0]):
#         smallestDistance = 0
#         for k in range(2):
#             distance = abs(sim_doublets[i] - centroids[k])
#             if k == 0:
#                 smallestDistance = distance
#                 classifications[i] = k
#             elif distance < smallestDistance:
#                 smallestDistance = distance
#                 classifications[i] = k

###Barcode and Doublet Prediction File
#https://github.com/swolock/scrublet/issues/5
outputfile=sample+"_scrublet_output_table.csv"
df = pd.DataFrame({
        'doublet_score': scrub.doublet_scores_obs_,
        'predicted_doublet': scrub.predicted_doublets_
})
barcodes.columns = ['Barcodes']
dffinal=pd.concat([barcodes,df],axis=1)
dffinal.to_csv(outputfile, index=False)

###DOUBLET HISTOGRAM and PREDICTED DOUBLET RATE
plt.figure(figsize=(3, 3))
scrub.plot_histogram()
plt.savefig(sample+"_doublets_hist.pdf")
plt.close()

print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
print('Done.')

plt.figure(figsize=(3, 3))
scrub.plot_embedding('UMAP', order_points=True)
plt.savefig(sample+"_doublets_umap.pdf")
plt.close()
