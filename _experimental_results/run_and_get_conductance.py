# Please follow the instructions in the README to install the python library interface and
# and copy mtkahypar.so to this folder to run the examples.
# use sudo (Linux & MacOS) or run shell as an administrator (Windows) to install system-wide
# make install-mtkahypar  

import os
import multiprocessing
import mtkahypar
import time

scaling_factor = 2147483647 / 1000
default_path = "/home/mmeshchaninova/mt-kahypar/iap_mt_kahypar/tests/interface/test_instances/ibm01.hgr"
# default_path = "/home/graph_collection/henrik/data/hypergraphs/talg_heuer/medium/ibm01.hgr"

def getNumBlocks():
  if input("Do you want to use a custom number of blocks? (y/n): ").lower() == 'y':
    num_blocks = int(input("Enter the number of blocks: "))
  else:
    num_blocks = 2
    print("Using default number of blocks: 2")
  return num_blocks

def getEpsilon():
  if input("Do you want to use a custom epsilon? (y/n): ").lower() == 'y':
    epsilon = float(input("Enter the value of epsilon: "))
  else:
    epsilon = 0.03
    print("Using default value of epsilon: 0.03")
  return epsilon

def getObjective():
  if input("Do you want to use a custom objective? (y/n): ").lower() == 'y':
    print("Available objectives: ")
    print("1. conductance_local")
    print("2. cut")
    print("3. km1")
    print("4. conductance_global")
    obj = str(input("Enter the number of the objective: "))
    if obj == "1":
      return mtkahypar.Objective.CONDUCTANCE_LOCAL
    elif obj == "2":
      return mtkahypar.Objective.CUT
    elif obj == "3":
      return mtkahypar.Objective.KM1
    elif obj == "4":
      return mtkahypar.Objective.CONDUCTANCE_GLOBAL
    else:
      print("Invalid input. Using default objective: conductance_local")
      return mtkahypar.Objective.CONDUCTANCE_LOCAL
  else:
    print("Using default objective: conductance_local")
    return mtkahypar.Objective.CONDUCTANCE_LOCAL

def getHGPath():
  if input("Do you want to use the default hypergraph (y/n): ").lower() == 'y':
    path_to_file = default_path
    print(f"Using default file {path_to_file}")  
  else:
    if input("Do you want to use your own hypergraph? (y/n): ").lower() == 'y':
      path_to_file = input("Enter the path to the hypergraph file: ")
      while not os.path.exists(path_to_file):
        print(f"File {path_to_file} does not exist.")
        if input("Try again? y/n: ").lower() == 'y':
          path_to_file = input("Enter the path to the hypergraph file: ")
        else:
          print("Using default file { path_to_file.split("/")[-1] } instead")
          path_to_file = default_path
          break
    else:
      print("Using hypergraph from the default directory")
      path_to_dir = "/home/graph_collection/henrik/data/hypergraphs/talg_heuer/"
      if input("Do you want to use a medium hypergraph? (y/n): ").lower() == 'y':
        path_to_dir += "medium/"
      else:
        print("Using large hypergraphs") 
        path_to_dir += "large/"
      file_name = input("Enter the name of the hypergraph file: ")
      path_to_file = path_to_dir + file_name
      while not os.path.exists(path_to_file) or file_name == "":
        print(f"File {path_to_file} does not exist.")
        if input("Try again? y/n: ").lower() == 'y':
          file_name = input("Enter the name of the hypergraph file: ")
          path_to_file = path_to_dir + file_name
        else:
          print("Using default file ibm01.hgr instead")
          path_to_file = default_path
          break
  return path_to_file

def getConductance():
  total_volume = 0
  part_volumes = [0 for _ in range(context.k)]
  cut_weights = [0 for _ in range(context.k)]
  for he in range(0, hypergraph.num_edges()):
    blocks = [False for _ in range(context.k)]
    num_blocks = 0
    for hn in hypergraph.pins(he):
      total_volume += hypergraph.edge_weight(he)
      block_id_hn = partitioned_hg.block_id(hn)
      part_volumes[block_id_hn] += hypergraph.edge_weight(he)
      if blocks[block_id_hn] == False:
        blocks[block_id_hn] = True
        num_blocks += 1
    if num_blocks > 1:
      for i in range(0, context.k):
        if blocks[i] == True:
          cut_weights[i] += hypergraph.edge_weight(he)
  max_conductance = -1
  num = -1
  denom = 1
  for i in range(0, context.k):
    min_vol = min(part_volumes[i], total_volume - part_volumes[i])
    if min_vol == 0:
      continue
    concuctance = cut_weights[i] / min_vol
    if max_conductance < concuctance:
      max_conductance = concuctance
      num = cut_weights[i]
      denom = min_vol
  return max_conductance, num, denom

def checkConductance(partitioned_hg):
  if input("Do you want to check the conductance? (y/n): ").lower() == 'y':
    max_conductance, num, denom = getConductance()
    print("actual conductance = " + str(max_conductance))
    print("conductance = " + str(partitioned_hg.conductance_local()))
    print("conductance = " + str(max_conductance * scaling_factor))
    if (num * scaling_factor / denom // 1) == partitioned_hg.conductance_local():
      print("\033[1;92mConductance was calculated correctly :)\033[0m")

mydir = os.path.dirname(os.path.realpath(__file__))

# Initialize
# mtk = mtkahypar.initialize(multiprocessing.cpu_count()) # use all available cores
mtk = mtkahypar.initialize(6) # use 6 cores

# Setup partitioning context
context = mtk.context_from_preset(mtkahypar.PresetType.CLUSTER) # use clustering preset
# In the following, we partition a hypergraph into four blocks
# with an allowed imbalance of 3% and optimize the connectivity metric

num_blocks = getNumBlocks()
objective = getObjective()
epsilon = getEpsilon()


context.set_partitioning_parameters(
  num_blocks,                       # number of blocks
  epsilon,                    # imbalance parameter
  objective ) # objective function
mtkahypar.set_seed(42)      # seed
context.logging = True

path_to_file = getHGPath()

# Load hypergraph from file
hypergraph = mtk.hypergraph_from_file(
  path_to_file, # hypergraph file
  context,
  mtkahypar.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Add fixed vertices
# hypergraph.add_fixed_vertices_from_file(
#   mydir + "/../tests/interface/test_instances/ibm01.k4.p1.fix", 4)

# Partition hypergraph
# start precise timer from time
start = time.time()
partitioned_hg = hypergraph.partition(context)
end = time.time()
print("Partitioning time = " + str(end- start) + " seconds")

# Output metrics
print("Partition Stats:")
print("Imbalance = " + str(partitioned_hg.imbalance(context)))
print("km1       = " + str(partitioned_hg.km1()))
print("cut       = " + str(partitioned_hg.cut()))
print("conductance = " + str(partitioned_hg.conductance_local() / scaling_factor))

checkConductance(partitioned_hg)