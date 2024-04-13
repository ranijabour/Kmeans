import argparse
from  enum import Enum
import pandas as pd
import numpy as np
import spkmeansmod

class Goal(Enum):
    SPK = 'spk'
    WAM = 'wam'
    DDG = 'ddg'
    LNORM = 'lnorm'
    JACOBI = 'jacobi'

def getFirstKCentroids(k,points):
    nrows=points.shape[0]
    arraylist = points.tolist()
    i=1
    np.random.seed(0)
    random_index= np.random.choice(nrows)
    centroids=[]
    centroid_indexes=[]
    centroids.append(arraylist[random_index])
    centroid_indexes.append(random_index)
    print(random_index, end="")
    print(",", end="")
    while(i!=k):
        Dl_list=[]
        P=[]
        sumDl=0
        for l in range(0,nrows):
            minDl=pow(np.linalg.norm(points[l]-centroids[0]),2)
            for j in range(1,i):
                Dl = pow(np.linalg.norm(points[l] - centroids[j]), 2)
                if (Dl < minDl):
                    minDl = Dl
            Dl_list.append(minDl)
            sumDl+=minDl
        for l in range(0,nrows):
            Px=Dl_list[l]/sumDl
            P.append(Px)
        i+=1
        m=np.random.choice(nrows, p=P)
        centroid_indexes.append(m)
        centroids.append(arraylist[m])
        print(m, end="")
        if (i< k):
            print(",", end="")
    print("")
    return centroids

def printCentroids(centroidsfromc):
    for i in range(k):
        output = ''
        for j in range(len(centroidsfromc[i])):
            if (j == 0):
                output = output +str(format(centroidsfromc[i][j],".4f"))
                output = output + ','
            else:
                if (j == (len(centroidsfromc[i]) - 1)):
                    output = output + str(format(centroidsfromc[i][j],".4f"))
                else:
                    output = output + str(format(centroidsfromc[i][j],".4f"))
                    output = output + ','
        if(i==k-1):
            print(output)
        else:
            print(output)

def print_mat(mat,rows,cols):
    for i in range(rows):
        output = ''
        for j in range(cols):
            if (j <cols-1):
                output = output + str(format(mat[i][j], ".4f"))
                output = output + ','
            else:
                output = output + str(format(mat[i][j], ".4f"))
        print(output)


# using ArgumentParser to initialize the parser so that we can add custom arguments
parser = argparse.ArgumentParser()
parser.add_argument("k",type=int)
parser.add_argument('goal', type=Goal, choices=list(Goal))
parser.add_argument("filename",type=str)
args = parser.parse_args()
k = args.k
assert k>=0, "Invalid Input!"
assert k!=1,"Invalid Input!"
goal = args.goal
filename = args.filename
df = pd.read_csv(filename,header=None)
points=df.to_numpy()
Points_aslist = points.tolist()
assert k<len(points)


if goal==Goal.WAM:
    res_mat=spkmeansmod.fitoperation( Points_aslist,k, len(Points_aslist), len(Points_aslist[0]),0)
    print_mat(np.asarray(res_mat),len(Points_aslist),len(Points_aslist))
elif goal == Goal.DDG:
    res_mat=spkmeansmod.fitoperation(Points_aslist,k, len(Points_aslist), len(Points_aslist[0]),1)
    print_mat(np.asarray(res_mat),len(Points_aslist),len(Points_aslist))
elif goal==Goal.LNORM:
    res_mat=spkmeansmod.fitoperation(Points_aslist,k, len(Points_aslist), len(Points_aslist[0]),2)
    print_mat(np.asarray(res_mat),len(Points_aslist),len(Points_aslist))
elif goal==Goal.JACOBI:
    res_mat=spkmeansmod.fitoperation(Points_aslist,k, len(Points_aslist), len(Points_aslist[0]),3)
    print_mat(np.asarray(res_mat),len(Points_aslist)+1,len(Points_aslist))
else:
    if goal == Goal.SPK:
        Tmatrix = spkmeansmod.fit_getTmat(Points_aslist, k, len(Points_aslist), len(Points_aslist[0]))
        Tmatrixasarray = np.array(Tmatrix)
        k = len(Tmatrixasarray[0])
        centroids = getFirstKCentroids(k, Tmatrixasarray)
        results_spk = spkmeansmod.fit_spk(Tmatrix, centroids, k, len(Points_aslist), k)
        printCentroids(np.array(results_spk))
