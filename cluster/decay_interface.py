import sys

sys.path.append("../")
import argparse

from grf.decay import DecaySpectra
from grf.pk_interp import PowerSpectrumGridInterpolator
from grf.units import *
from tqdm import *

parser = argparse.ArgumentParser()

parser.add_argument("--i_m_Ap", action="store", dest="i_m_Ap", default=0, type=int)
parser.add_argument("--use_stellar", action="store", dest="use_stellar", default=0, type=int)
parser.add_argument("--use_arcade2", action="store", dest="use_arcade2", default=0, type=int)
parser.add_argument("--use_firas", action="store", dest="use_firas", default=0, type=int)

results = parser.parse_args()

i_m_Ap = results.i_m_Ap
use_stellar = results.use_stellar
use_arcade2 = results.use_arcade2
use_firas = results.use_firas

tag = str(i_m_Ap) + "_" + str(use_stellar) + str(use_arcade2) + str(use_firas)

output_dir = "/scratch/sm8383/chi2_arys/decay_"

pspec = PowerSpectrumGridInterpolator("franken_lower")
spec_scan = DecaySpectra(pspec)

m_Ap = 1e-12 * eV
m_a = 1e-6 * eV

tau_a = spec_scan.tau_a_LL()

m_a_ary = np.logspace(-6, -1, 50) * eV
m_Ap_ary = np.logspace(-14, -9, 50) * eV

m_Ap = m_Ap_ary[i_m_Ap]

ratio_ary = []

if use_firas:
    for m_a in tqdm(m_a_ary):
        z_ary, P_ary, max_rat, eps_max = spec_scan.get_max_CMB_photon_ratio_FIRAS(m_Ap, m_a, tau_a, z_at_which_max=17, one_plus_delta_bound=1e2)
        ratio_ary.append(max_rat)
else:
    for m_a in tqdm(m_a_ary):
        z_ary, P_ary, max_rat, eps_max = spec_scan.get_max_CMB_photon_ratio_z_17(m_Ap, m_a, tau_a, z_at_which_max=17, one_plus_delta_bound=1e2, use_stellar=use_stellar, use_arcade2=use_arcade2)
        ratio_ary.append(max_rat)

np.savez(output_dir + tag, ratio_ary=ratio_ary)
