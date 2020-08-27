import sys

sys.path.append("../")

import os
import random
import numpy as np

from grf.units import *

batch = """#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH -t 10:00:00
#SBATCH --mem=4GB
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=sm8383@nyu.edu

cd /home/sm8383/spectral-distortions-perturbations/cluster

conda activate

"""

# Main set of plots

i_m_Ap_ary = np.arange(50)

for i_m_Ap in i_m_Ap_ary:
    batchn = batch + "\n"
    batchn += "python decay_interface.py --i_m_Ap " + str(i_m_Ap) + " --use_firas 1"
    fname = "batch/submit.batch"
    f = open(fname, "w")
    f.write(batchn)
    f.close()
    os.system("chmod +x " + fname)
    os.system("sbatch " + fname)

# for i_m_Ap in i_m_Ap_ary:
#     batchn = batch + "\n"I
#     batchn += "python decay_interface.py --i_m_Ap " + str(i_m_Ap) + " --use_arcade2 1"
#     fname = "batch/submit.batch"
#     f = open(fname, "w")
#     f.write(batchn)
#     f.close()
#     os.system("chmod +x " + fname)
#     os.system("sbatch " + fname)

# for i_m_Ap in i_m_Ap_ary:
#     batchn = batch + "\n"
#     batchn += "python decay_interface.py --i_m_Ap " + str(i_m_Ap) + " --use_stellar 1"
#     fname = "batch/submit.batch"
#     f = open(fname, "w")
#     f.write(batchn)
#     f.close()
#     os.system("chmod +x " + fname)
#     os.system("sbatch " + fname)
