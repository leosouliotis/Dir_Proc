# Dir_Proc
Dirichlet process collapsed Gibbs for infinite Gaussian mixtures

# Dependencies and prerequisites
* Julia Lang (v>0.6) with the following packages installed

  * DataFrames
  * Distributions
  * StatsFuns
  * Clustering
  * ArgParse

More information on how to install Julia in your machine can be found [here](https://julialang.org/downloads/).
The prerequisite packages can be downloaded the following way:

  ```julia
   Pkg.add("DataFrames")
   Pkg.add("Distributions")
   Pkg.add("StatsFuns")
   Pkg.add("Clustering")
   Pkg.add("ArgParse")
   ```

# Infinite Gaussian mixtures

The prior mean vector **μ0** and the prior covariance matrix **Λ0** have been set to a d dimensional **0** vector and a diagonal matrix with entries 1/(d * w), where w is the diagonal elements of the covariance matrix of the data, where d is the dimension of the data set.
*User is free to change these values in the script*.

 The parameter α that control the assignment of the data in a new cluster is set to 1. More information about DP can be found [here](https://www.stats.ox.ac.uk/~teh/research/npbayes/Teh2010a.pdf).
 
 # Instructions
 
 The DP clustering can be performed as

`julia DP_CS.jl -d Dataset -i No_of_iterations -k Initial_number_of_clusters`

where
  * -d : Dataset of interest (N x d)
  * -i : Number of iteration for the Gibbs sampler
  * -k : Initial number of clusters

The results of the clustering are saved in a .csv file named DP_clustering.csv

# Example

We simulated a 5 dimensional dataset consisting of 4 clusters. Each cluster has a different mean vector and the same covariance matrix. After downloading the data and the true cluster allocation, you can run the DP Gibbs Sampler as

``julia DP_CS.jl -d Sample_data.csv -i 15 -k 10``

starting with 10 clusters for 15 iterations.

In order to validate our clustering, knowing the true cluster allocation of the observations, we use the following script which calculates the Rand Index between the true cluster allocation and the allocation under the DP model

`using DataFrame
using Clustering

rca = readtable("sample_rca.csv",header= false)
dpca = readtable("DP_clustering.csv",header= false)

randindex(Array(rca),Array(dpca))[1]
`
The Rand Index is equal to 1, which results that the clustering under the DP model is perfect
