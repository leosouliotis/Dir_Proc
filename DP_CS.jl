#!/usr/bin/env julia

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
    "--opt2", "-k"
      arg_type = Int
      help = "Starting number of clusters"
    "--opt3", "-i"
      arg_type = Int
      help = "No of Iterations"
  end
  return parse_args(s)
end

# Import and sort coverage dataset
parsed_args = parse_commandline()
data = readtable(parsed_args["opt1"],header = true)
N,d = size(data)
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
K = Int(parsed_args["opt2"])
label = kmeans(transpose(X),K,maxiter=1000, tol=0.00001).assignments
println("Initial ARI:")
println(randindex(rga,label)[1])
nk = counts(label);
prob = zeros(Float64,K+1);
println("Gibbs sampler initialized and  starting...")

# Run Gibbs Sampler
Iter = Int(parsed_args["opt3"])
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
end
writedlm("DP_clustering_c10K_test.csv", label, ",")
