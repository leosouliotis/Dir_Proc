using DataFrames
using Distributions
using StatsFuns
using Clustering
using ArgParse
###################################################
function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--opt1" ,"-d"
      help = "A CSV file containing youd data"
    "--opt2", "-o"
      help = "True Cluster Allocation"
    "--opt3", "-k"
      arg_type = Int
      help = "Starting number of clusters"
    "--opt4", "-i"
      arg_type = Int
      help = "No of Iterations"
  end
  return parse_args(s)
end

# Import and sort coverage dataset
parsed_args = parse_commandline()
data = readtable(parsed_args["opt1"],header = true)
data = sort(data,cols=:contig_id)
contigs= data[:,2]
data = delete!(data, :1)
data = delete!(data, :1)
N,d = size(data)

# Import and sort trur cluster allocation
tga=readtable(parsed_args["opt2"],header = true)
rga = sort(tga,cols=:clustering_species_1_)
println(rga)
rga=tga[:,3]
println("Data imported")


# Multivariate t distribution density calulation
function dmvtp(x, mu, Sigma, df)
  p= size(X)[2]
  Omega = inv(Sigma)
  ss = x - vec(mu)
  z = ss'*Omega*ss
  dens = lgamma((df + p)/2) - (lgamma(df/2) + (p/2) *
    log(df) + (p/2) * log(π) + 1/2 * logdet(Sigma) + (df +
    p)/2 * log(1 + (1/df) * z))
end

# Hyperparameters
X = Matrix(data)
Λ_0 = inv(diagm(1./(d*diag(cov(X)))))
μ_0 = zeros(Float64,(d))
κ_0 = 0.001
ν_0 = d
α = 1

# Initialization
K = Int(parsed_args["opt3"])
label = kmeans(transpose(X),K,maxiter=1000, tol=0.00001).assignments
println("Initial ARI:")
println(randindex(rga,label)[1])
nk = counts(label);
prob = zeros(Float64,K+1);
println("Gibbs sampler initialized and  starting...")

# Run Gibbs Sampler
rand_prev = rand_diff = 1
#iter = 1
#while abs(rand_diff) > 0.01
for iter in 1:30
  for datum = 1:N
    x_i = X[datum,:]
    z_i = label[datum]
    label[datum] = 0
    nk[z_i] -= 1
    if in(z_i,label) == false
       K -= 1
       deleteat!(nk,z_i)
       ind = find(label.>z_i)
       label[ind] -= 1
    end
    prob = zeros(Float64,K+1)
    for k = 1:K
      ind = find(label.==k)
      data_sub = X[ind,:]
      κ_n = κ_0 + nk[k] - 1
      ν_n = ν_0 + nk[k] - 1
      mc = mean(data_sub,1)
      μ = (κ_0*μ_0 + nk[k]*vec(mc))/κ_n
      Λ_n = Λ_0 + scattermat(data_sub) + (κ_0*nk[k])/κ_n * (vec(mc)-μ_0)*(vec(mc)-μ_0)'
      prob[k] = log(nk[k]) + dmvtp(x_i, μ, (κ_n + 1)/(κ_n*(ν_n - d +1))*Λ_n, ν_n - d +1)[1]
    end
    prob[K+1] = log(1) + dmvtp(x_i , μ_0, (κ_0 + 1)/(κ_0*(ν_0 - d +1))*Λ_0, ν_0 - d +1)[1]
    k = rand(Categorical(exp.(prob-logsumexp(prob))))
    if k == K+1
        append!(nk,1)
        K += 1
    else
        nk[k] += 1
    end
    label[datum] = k
  end
  #rand_cur = randindex(rga,label)[1]
  #println("After Iter $iter, the ARI is $rand_cur with $K clusters")
  println("After Iter $iter, $K clusters")
  #iter = iter + 1
  #rand_dif = rand_cur - rand_prev
  #rand_prev = rand_cur
  #writedlm("DP_clustering_c10K_test.csv", label, ",")
end

'''
# for all clusters estimated by the algorithm
for ind in unique(label)
  println("Observed cluster: $ind")
  temp = find(label.==ind)
  println("Contigs in observer cluster $ind: ", length(temp))
  println("--------------------------------")
  # find the number of estimated clusters in the original clusters
  t = unique(rga[temp])
  temp_cum = 0
  for j in t
    temp_in = length(find(rga[temp].==j))
    println("True cluster $j ","|",temp_in)
    temp_cum += temp_in
  end
  println("\n")
end
'''
label = readtable("DP_clustering_c10K_test.csv",header = false)
label = label[:,1]
# Find 'pure' and 'dirty clusters'
for ind in unique(rga)
  println("True cluster: $ind")
  temp = find(rga.==ind)
  println("Contigs in true cluster $ind: ", length(temp))
  println("--------------------------------")
  # find the number of estimated clusters in the original clusters
  t = unique(label[temp])
  temp_cum = 0
  for j in t
    temp_in = length(find(label[temp].==j))
    println("Observed cluster $j ","|",temp_in)
    temp_cum += temp_in
  end
  println("\n")
end


cl_graph = loadgraphs("Graph_fin_extra_links_test.xml.gz",GraphMLFormat())
cl_graph = cl_graph["G"]

t = length(find(rga.==17))
temp = zeros(Int64,(t,2))
ind = 1
for v in vertices(cl_graph)
  k = rga[v]
  if k == 17
    temp[ind,1] = length(out_neighbors(cl_graph,v))
    ind += 1
  end
end
print(counts(temp))
