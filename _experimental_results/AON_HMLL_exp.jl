###############################################################
# Script to run AON HMLL on simulated data and analyze results
# Source: https://github.com/veronicapoda/modularity/tree/main
# Original file name: AON_HMLL_script.jl
###############################################################


###### Install packages if needed
#using Pkg
#Pkg.add("HyperModularity")
#Pkg.add("SimpleHypergraphs")

using HyperModularity
using SimpleHypergraphs 
using DelimitedFiles
using ArgParse
using Clustering

function parse_commandline()
    s = ArgParse.ArgParseSettings()

    ArgParse.@add_arg_table! s begin
        "-d"
            arg_type = String
            default = "./"
            help = "data directory"
        "-m"
            arg_type = String
            required = true
            help = "generating model"
        "-s"
            arg_type = String
            required = true
            help = "scenario (parameter choice)"
        "-n"
            arg_type = Int64
            default = 25
            help  = "num of hypergraphs in the scenario"
        "-i"
            arg_type = String
            default = ""
            help  = "which instances to use: e.g. 1,2,5"
        "-v"
            arg_type = String
            default = "AON_HMLL_clique"
            help = "version of mt-kahypar used (for output file naming)"
        
    end
    s.usage = "julia AON_HMLL_script.jl -d <data_root_dir> -m <model> -s <scen> -i \"1,2,3\""
    return ArgParse.parse_args(s)
end
######

# exemple of a command
# julia AON_HMLL_script.jl -v "test" -d "survey_benchmark/" -m "HyperSBM/" -s "scenA1/" -i "1,2,3" -t 3


############
### Quality of clustering
############

# Adjusted Rand Index
function ari(x,y)
    evaluations = randindex(x, y)
    ari = evaluations[1]
    return ari
end

############
### Function from Chodrow et al. to read data
### (slightly modified to accept more general forms of dataname)
############
function read_hypergraph_edges(dataname::String, maxsize::Int64=1000, minsize::Int64=2)
    E = Dict{Integer, Dict}()
    open(string(dataname,"_he.txt")) do f
        for line in eachline(f)
            edge = [parse(Int64, v) for v in split(line, ',')]
            sort!(edge)
            if minsize <= length(edge) <= maxsize
                sz = length(edge)
                if !haskey(E, sz)
                    E[sz] = Dict{}()
                end
                E[sz][edge] = 1
            end
        end
    end
    return E
end

function writePartition(filename::String, partition::Vector{Int})
    open(filename, "w") do f
        for v in partition
            write(f, string(v)*"\n")
        end
    end
end


#######################
#######################
function main()
    parsed_args = parse_commandline()
    
    ###### PARSE ARGUMENTS ######
    data_root_dir = parsed_args["d"]
    cd(data_root_dir)
    model = parsed_args["m"]
    scenario = parsed_args["s"]
    nreps =  parsed_args["n"]
    version = parsed_args["v"]
    instances_str = parsed_args["i"]
    if instances_str != ""
        instances = parse.(Int64, split(instances_str, ","))
        nreps = length(instances)
    else 
        instances = 1:nreps
    end

    
    ###### INITIALIZATIONS OF VARIABLES ######
    Ari = zeros(Float64, nreps)
    runtimes = zeros(Float64, nreps)
    Q = zeros(Float64, nreps)
    Q_true = zeros(Float64, nreps)
    K_true = zeros(Int64, nreps)
    K_hat = zeros(Int64, nreps)
    relative_K_error = zeros(Float64, nreps)
    relative_mod_error = zeros(Float64, nreps)

    ###### LOOP OVER DATASETS ######
    for i in 1:nreps
        instance = instances[i]
        # print("\n Instance: $i \n")
                
        # read true label vector
        filename = string(model,scenario,"rep",instance,"_assign.txt")
        Z_true = vec(readdlm(filename, Int))
        K_true[i] = length(unique(Z_true)) # true number of clusters
        
        runtimes[i] = @elapsed begin
            # read data
            filename = string(model,scenario,"rep",instance)
            E = read_hypergraph_edges(filename) # possible to add min and max hyperedge size
            n = maximum([maximum(e) for k in keys(E) for e in keys(E[k])])
            D = zeros(Int64, n)
            for (sz, edges) in E
                for (e, _) in edges
                    D[e] .+= 1
                end
            end
            maxedges = maximum(keys(E))
            for k in 1:maxedges
                if !haskey(E, k)
                    E[k] = Dict{}()
                end
            end
            # construct hypergraph
            N = 1:n
            H = hypergraph(N, E, D)
            
            ############ AON HMLL ################
            Z_hat = Simple_AON_Louvain(H, startclusters = "cliquelouvain") # By default randflag::Bool=false - the result is non random
        end

        # performance measures
        K_hat[i]= length(unique(Z_hat)) # estimated number of clusters
        relative_K_error[i] = (K_true[i] - K_hat[i]) / K_true[i]

        ######## Extract omega parameters to compute modularity ######
        # compute initialisation corresponding to startclusters = "cliquelouvain"
        Z0 = HyperModularity.CliqueExpansionModularity(H,1.0) # default value gamma=1.0
        He2n, weights = HyperModularity.hypergraph2incidence(H)
        e2n = incidence2elist(He2n)
        d = 1.0*H.D # degree vector - coercion to float for compatibility with AON_louvain
        β, γ,omega = learn_omega_aon(e2n,weights,Z0,maxedges,d,n)

        Q[i] = modularity_aon(H,Z_hat,omega)
        Q_true[i] = modularity_aon(H,Z_true,omega)
        relative_mod_error[i] = (Q_true[i] - Q[i])/(Q_true[i])

        writePartition(string(model,scenario,version,"/rep",instance,"_he.hgr.part"), Z_hat)
       
        ######### Compute ARI #######
        Ari[i] = ari(Z_hat,Z_true)
        
        # save partial results
        filename= string(model,scenario,"$(version).partial_results")
        open(filename, "w") do f
            write(f, "CPU time = $runtimes"*"\n"*"Modularity = $Q"*"\n"*"GT_Modularity = $Q_true"*"\n"*"Relative Modularity error = $relative_mod_error"*"\n"*"K = $K_hat"*"\n"*"GT_K = $K_true"*"\n"*"Relative k error = $relative_K_error"*"\n"*"Instances = $(instances)"*"\n")
            # write(f, "ARI = $Ari"*"\n"*"CPU time = $runtimes"*"\n"*"Modularity = $Q"*"\n"*"GT_Mod = $Q_true"*"\n"*"Relative Mod. error = $relative_mod_error"*"\n"*"K_hat = $K_hat"*"\n"*"K_true = $K_true"*"\n"*"Relative k error = $relative_K_error"*"\n"*"Instances = $(instances)")
            # write(f, "ARI=$(string(Ari))"*"\n"*"CPU time = $runtimes"*"\n"*"Modularity = $Q"*"\n"*"GT_Mod = $Q_true"*"\n"*"K_hat = $K_hat"*"\n"*"K_true = $K_true")
        end
    end
    
    ##########
    # Correct the first computing time
    instance = instances[1]
    runtimes[1] = @elapsed begin
        # read data
        filename = string(model,scenario,"rep",instance)
        E = read_hypergraph_edges(filename) # possible to add min and max hyperedge size
        n = maximum([maximum(e) for k in keys(E) for e in keys(E[k])])
        D = zeros(Int64, n)
        for (sz, edges) in E
            for (e, _) in edges
                D[e] .+= 1
            end
        end
        maxedges = maximum(keys(E))
        for k in 1:maxedges
            if !haskey(E, k)
                E[k] = Dict{}()
            end
        end
        # construct hypergraph
        N = 1:n
        H = hypergraph(N, E, D)
        
        ############ AON HMLL ################
        Z_hat = Simple_AON_Louvain(H, startclusters = "cliquelouvain") # By default randflag::Bool=false - the result is non random
    end
    #########
    

    ###### OUTPUT FILES #######
    filename= string(model,scenario,"$(version).results")
    open(filename, "w") do f
        write(f, "CPU time = $runtimes"*"\n"*"Modularity = $Q"*"\n"*"GT_Modularity = $Q_true"*"\n"*"Relative Modularity error = $relative_mod_error"*"\n"*"K = $K_hat"*"\n"*"GT_K = $K_true"*"\n"*"Relative k error = $relative_K_error"*"\n"*"Instances = $(instances)"*"\n")
        # write(f, "ARI = $Ari"*"\n"*"CPU time = $runtimes"*"\n"*"Modularity = $Q"*"\n"*"GT_Mod = $Q_true"*"\n"*"Relative Mod. error = $relative_mod_error"*"\n"*"K_hat = $K_hat"*"\n"*"K_true = $K_true"*"\n"*"Relative k error = $relative_K_error"*"\n"*"Instances = $(instances)")
    end

end


##################
##################
main()