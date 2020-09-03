import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.misc import derivative
import astropy.units as u
import astropy.constants as c
from astropy.cosmology import FlatLambdaCDM, z_at_value
from tqdm import *

from grf.grf import FIRAS
from grf.units import *


class DecaySpectra(FIRAS):
    def __init__(self, log_pk_b_interp_fn, cosmo=None, z_reio=7.82):
        """
        Class to compute CMB and 21cm spectra.
        :param cosmo: Cosmology as an astropy object
        :param z_reio: Redshift of reionization
        """
        FIRAS.__init__(self, log_pk_b_interp_fn, cosmo=cosmo, z_reio=z_reio)

        self.omega_21 = 5.9e-6 * eV  # Angular frequency of 21cm photons
        self.tau_u = self.cosmo.age(0).value * 1e9 * Year  # Age of the Universe
        self.rho_c = self.cosmo.critical_density0.value * Gram / Centimeter ** 3

        df_all = pd.read_csv("../data/apjlaabf86t2_ascii.txt", delimiter="\t", skiprows=list(range(6)) + list(range(21, 24)), names=["freq", "survey", "fid", "unc", "n/a"])
        df_arcade2 = df_all[df_all["survey"].str.contains("ARCADE")]
        df_llfss = df_all[df_all["survey"].str.contains("LLFSS")]
        df_radio = df_all[~df_all["survey"].str.contains("LLFSS|ARCADE")]

        self.omega_ARCADE2_ary = 2 * np.pi * df_arcade2["freq"].values * 1e9 / Sec
        self.T_ARCADE2_fid = df_arcade2["fid"].values
        self.T_ARCADE2_unc = df_arcade2["unc"].values

        self.omega_radio_ary = 2 * np.pi * df_radio["freq"].values * 1e9 / Sec
        self.T_radio_fid = df_radio["fid"].values
        self.T_radio_unc = df_radio["unc"].values

        self.omega_llfss_ary = 2 * np.pi * df_llfss["freq"].values * 1e9 / Sec
        self.T_llfss_fid = df_llfss["fid"].values
        self.T_llfss_unc = df_llfss["unc"].values

        self.omega_all_radio_ary = np.concatenate((self.omega_ARCADE2_ary, self.omega_radio_ary, self.omega_llfss_ary))

    def T_CMB_powerlaw(self, A_r, beta, T_CMB_0, omega_obs, omega_base=2 * np.pi * 78 * 1e6 / Sec / eV):
        """ Powerlaw spectrum according to Eq. (7) of 1902.02438
        """
        return T_CMB_0 * (1 + A_r * (omega_obs / omega_base) ** beta)

    def get_powerlaw_params(self):
        omega_ary = np.concatenate((self.omega_ARCADE2_ary / eV, self.omega_radio_ary / eV, self.omega_llfss_ary / eV))
        T_ary = np.concatenate((self.T_ARCADE2_fid, self.T_radio_fid, self.T_llfss_fid))
        T_unc_ary = np.concatenate((self.T_ARCADE2_unc, self.T_radio_unc, self.T_llfss_unc))

        def chi2(x):
            A_r, beta, T_CMB_0 = x
            return np.sum(((self.T_CMB_powerlaw(A_r, beta, T_CMB_0, omega_ary) - T_ary) / T_unc_ary) ** 2)

        self.opt_powerlaw = minimize(chi2, x0=[100, -2.6, 2.725])

    def T_CMB_powerlaw_ary(self, z_ary, f_excess=0.007):
        A_r_fit, beta_fit, T_CMB_0_fit = self.opt_powerlaw.x
        return self.T_CMB_powerlaw(f_excess * A_r_fit, beta_fit, T_CMB_0_fit, self.omega_21 / eV / (1 + z_ary)) * (1 + z_ary)

    def ratio_T_CMB_powerlaw(self, z, f_excess):
        A_r_fit, beta_fit, T_CMB_0_fit = self.opt_powerlaw.x
        return self.T_CMB_powerlaw(f_excess * A_r_fit, beta_fit, T_CMB_0_fit, self.omega_21 / eV / (1 + z)) * (1 + z) / (T_CMB_0_fit * (1 + z))

    def compute_temperature_evolution(self, m_a, tau_a_div_tau_u, z_ary, z_P_bin_centers, P_ary):
        """
        Computer 21cm temperature evolution, in standard scenario and with dark photon energy injection
        :param m_a: Mediator mass
        :param tau_a_div_tau_u: Mediator lifetime compared to lifetime of Universe
        :param z_ary: Redshifts at which to compute temperatures
        :param z_P_bins_ary: Redshift bins at which probabilities are specified
        :param P_ary: Probabilities corresponding to redshift bins `z_P_bins_ary`, times omega_0
        """

        # Lifetime
        self.tau_a = self.tau_u * tau_a_div_tau_u

        # Redshift bins at which transition probabilities specified
        # z_P_bin_centers = 10 ** ((np.log10(z_P_bins_ary[1:]) + np.log10(z_P_bins_ary[:-1])) / 2.0)

        # Differential number density (spectrum) of A' and CMB photons
        self.dn_SM_domega_ary = self.dn_CMB_domega(omega=self.omega_21, z=z_ary)
        omega_0_ary = self.omega_21 / (1 + z_ary)
        self.dn_A_domega_ary = np.array([self.dn_A_domega(m_a=m_a, z=z, z_res_ary=z_P_bin_centers, P_ary=P_ary / omega_0_ary[i_z], omega=self.omega_21, tau_a=self.tau_a) for i_z, z in enumerate(z_ary)])

        # Final temperature evolution
        self.T_CMB_SM = self.T_CMB_K(z_ary)  # Standard scenario
        self.T_CMB_A = self.T_CMB_K(z_ary) * ((self.dn_SM_domega_ary + self.dn_A_domega_ary) / self.dn_SM_domega_ary)  # With A' decays

    def T_CMB(self, z, T_0=None):
        """
        :param z: Redshift
        :return: CMB temperature in natural units
        """
        if T_0 is None:
            T_0_n = self.T_0_n
        else:
            T_0_n = (c.k_B * T_0 * u.Kelvin).to(u.eV).value * eV
        return T_0_n * (1 + z)

    def T_CMB_K(self, z):
        """
        :param z: Redshift
        :return: CMB temperature in K
        """
        return self.T_0 * (1 + z)

    def dn_A_domega(self, m_a, z, z_res_ary, P_ary, omega, tau_a, return_Ap_spec=False):
        """
        :param m_a: Mass of DM
        :param z: Redshift at which spectrum of A is calculated
        :param z_res_ary: Array or redshifts over which transition probability specified
        :param P_ary: Transition probabilities in `z_res_ary`
        :param omega: Frequency at z
        :param tau_a: DM lifetime
        :return: Differential number density of photon A, dn/omega at frequency omega
        """
        omega_0 = omega / (1 + z)

        self.z_dec = m_a / 2 / omega_0 - 1
        self.t_c = self.cosmo.age(self.z_dec).value * 1e9 * Year
        self.H = self.cosmo.H(self.z_dec).value * Kmps / Mpc

        omega_res_ary = omega_0 * (1 + z_res_ary)

        dn_Ap_domega_0 = 2 * self.cosmo.Odm0 * self.rho_c / (m_a * tau_a * omega_0 * self.H) * np.heaviside(self.z_dec, 0) * np.exp(-self.t_c / tau_a)

        P_tot = np.sum(P_ary * np.heaviside(z_res_ary - z, 0) * np.heaviside(m_a / 2 / omega_res_ary - 1, 0))
        dn_Ap_domega = (1 + z) ** 2 * dn_Ap_domega_0

        dn_A_domega = dn_Ap_domega * P_tot

        if return_Ap_spec:
            return dn_A_domega, dn_Ap_domega
        else:
            return dn_A_domega

    def dn_CMB_domega(self, omega, z=0, T_0=None):
        """
        :param omega: Photon frequency
        :param z: Redshift
        :return: Number density spectrum of CMB
        """
        return 8 * np.pi * (omega / (2 * np.pi)) ** 2 * 1 / (np.exp(omega / self.T_CMB(z, T_0)) - 1) / (2 * np.pi)

    def f_a(self, tau_a, m_a):
        """
        From Eq. 5 of 1803.07048
        :param tau_a: Lifetime
        :param m_a: Dark matter mass
        :return: f_a
        """
        return np.sqrt(m_a ** 3 / (64 * np.pi * tau_a ** -1))

    def eps_stellar_cooling(self, f_a):
        """ Maximum allowed coupling from stellar cooling bounds, from 1803.07048
        """
        return 2e-9 * GeV ** -1 * f_a

    def tau_a_LL(self):
        """ Lower limit on mediator lifetime, from 1606.02073
        """
        return 1 / ((6.3e-3 / 1e9) * Year ** -1)

    def chi2_FIRAS_extra_B(self, x, eps, P_tot, dn_A_domega_ary):
        """ The FIRAS \chi^2 statistic, to do minimization over blackbody temperature
        """
        T_0 = x[0]

        # t = (self.B_CMB(self.omega_FIRAS, T_0) * (1 - P_tot * (eps / self.eps_base) ** 2) + dn_A_domega_ary * (eps / self.eps_base) ** 2 * (self.omega_FIRAS / (4 * np.pi)) * (2 * np.pi)) / (1e6 * Jy)
        t = (self.B_CMB(self.omega_FIRAS, T_0) * (1 - P_tot * (eps / self.eps_base) ** 2)) / (1e6 * Jy)

        return np.dot((self.d - t), np.matmul(self.Cinv, (self.d - t)))

    def get_max_CMB_photon_ratio_FIRAS(self, m_Ap, m_a, tau_a, z_at_which_max=17, **kwargs):

        z_ary, dP_dz_ary, _, _, z_ary_homo, P_ary_homo = self.P_tot_perturb(self.omega_FIRAS, self.eps_base, m_Ap, get_total_prob=1, Ap_DM=0, one_plus_delta_bound=1e2, z_excise_min=6, z_excise_max=20)

        P_ary = (dP_dz_ary[:, 1:] + dP_dz_ary[:, :-1]) / 2.0 * np.diff(z_ary)
        z_res_ary = (z_ary[1:] + z_ary[:-1]) / 2.0

        dn_A_domega_ary = []

        for i_omega, omega in enumerate(self.omega_FIRAS):
            dn_A_domega = self.dn_A_domega(m_a=m_a, z=0, z_res_ary=np.concatenate([z_res_ary, z_ary_homo[i_omega]]), P_ary=np.concatenate([P_ary[i_omega], P_ary_homo[i_omega]]), omega=omega, tau_a=tau_a)
            dn_A_domega_ary.append(dn_A_domega)

        dn_A_domega_ary = np.array(dn_A_domega_ary)
        P_tot = np.array([np.sum(P_ary[i_omega]) + np.sum(P_ary_homo[i_omega]) for i_omega in range(len(self.omega_FIRAS))])

        self.eps_ary = np.logspace(-11.0, -1.0, 100)

        chi2_null = minimize(self.chi2_FIRAS, x0=[2.725], args=(0, np.ones_like(self.omega_FIRAS)), method="SLSQP").fun
        chi2_ary = np.zeros(len(self.eps_ary))

        for i_eps, eps in enumerate((self.eps_ary)):
            chi2_min = minimize(self.chi2_FIRAS_extra_B, x0=[2.725], args=(eps, P_tot, dn_A_domega_ary), method="SLSQP")
            chi2_ary[i_eps] = chi2_min.fun

        eps_max = self.get_lim([2 * (chi2_ary - chi2_null)])

        dn_CMB_domega = self.dn_CMB_domega(omega=self.omega_21, z=z_at_which_max)

        omega_0 = self.omega_21 / (1 + z_at_which_max)

        z_ary, dP_dz_ary, _, _, z_ary_homo, P_ary_homo = self.P_tot_perturb([omega_0], eps_max, m_Ap, z_excise_min=6, z_excise_max=20, get_total_prob=1, Ap_DM=0, **kwargs)
        dP_dz_ary = np.nan_to_num(dP_dz_ary[0])

        P_ary = (dP_dz_ary[1:] + dP_dz_ary[:-1]) / 2.0 * np.diff(z_ary)
        z_res_ary = (z_ary[1:] + z_ary[:-1]) / 2.0

        dn_A_domega = self.dn_A_domega(m_a=m_a, z=z_at_which_max, z_res_ary=np.concatenate([z_res_ary, z_ary_homo[0]]), P_ary=np.concatenate([P_ary, P_ary_homo[0]]), omega=self.omega_21, tau_a=tau_a)

        return np.concatenate([z_res_ary, z_ary_homo[0]]), np.concatenate([P_ary, P_ary_homo[0]]), ((dn_A_domega + dn_CMB_domega) / dn_CMB_domega), eps_max

    def get_max_CMB_photon_ratio_z_17(self, m_Ap, m_a, tau_a, use_stellar=1, use_arcade2=1, z_at_which_max=17, sigma=2.0, **kwargs):

        eps_base = 1e-6

        eps_max_list = []

        if use_arcade2:

            # List of frequencies of radio surveys
            omega_check_ary = np.concatenate([self.omega_radio_ary, self.omega_llfss_ary, [self.omega_ARCADE2_ary[0]]])

            # Upper 2\sigma temperatures from radio surveys
            T0_check_ary = np.concatenate([self.T_radio_fid + sigma * self.T_radio_unc, self.T_llfss_fid + sigma * self.T_llfss_unc, [(self.T_ARCADE2_fid + sigma * self.T_ARCADE2_unc)[0]]])

            # Corresponding upper limit on photon flux from radio surveys
            dn_CMB_domega_ARCADE2_upper_check = self.dn_CMB_domega(omega_check_ary, z=0, T_0=T0_check_ary) / (eV ** -1)

            dn_A_domega_ary = []

            z_ary, dP_dz_ary, _, _, z_ary_homo, P_ary_homo = self.P_tot_perturb(omega_check_ary, eps_base, m_Ap, get_total_prob=1, Ap_DM=0, z_excise_min=6, z_excise_max=20, **kwargs)

            P_ary = (dP_dz_ary[:, 1:] + dP_dz_ary[:, :-1]) / 2.0 * np.diff(z_ary)
            z_res_ary = (z_ary[1:] + z_ary[:-1]) / 2.0

            for i_omega, omega in enumerate(omega_check_ary):
                dn_A_domega = self.dn_A_domega(m_a=m_a, z=0, z_res_ary=np.concatenate([z_res_ary, z_ary_homo[i_omega]]), P_ary=np.concatenate([P_ary[i_omega], P_ary_homo[i_omega]]), omega=omega, tau_a=tau_a) / (eV ** -1)
                dn_A_domega_ary.append(dn_A_domega)

            dn_A_domega_ary = np.array(dn_A_domega_ary)

            eps_max_ARCADE2 = eps_base * np.sqrt(np.min(dn_CMB_domega_ARCADE2_upper_check / dn_A_domega_ary))

            eps_max_list.append(eps_max_ARCADE2)

        if use_stellar:
            f_a = self.f_a(tau_a, m_a)
            eps_max_stellar = self.eps_stellar_cooling(f_a)
            eps_max_list.append(eps_max_stellar)

        eps_max = np.min(eps_max_list)

        # Get ratio of 21-cm photon number densities at z=17

        dn_CMB_domega = self.dn_CMB_domega(omega=self.omega_21, z=z_at_which_max)

        omega_0 = self.omega_21 / (1 + z_at_which_max)

        z_ary, dP_dz_ary, _, _, z_ary_homo, P_ary_homo = self.P_tot_perturb([omega_0], eps_max, m_Ap, get_total_prob=1, Ap_DM=0, z_excise_min=6, z_excise_max=20, **kwargs)
        dP_dz_ary = dP_dz_ary[0]

        P_ary = (dP_dz_ary[1:] + dP_dz_ary[:-1]) / 2.0 * np.diff(z_ary)
        z_res_ary = (z_ary[1:] + z_ary[:-1]) / 2.0

        dn_A_domega = self.dn_A_domega(m_a=m_a, z=z_at_which_max, z_res_ary=np.concatenate([z_res_ary, z_ary_homo[0]]), P_ary=np.concatenate([P_ary, P_ary_homo[0]]), omega=self.omega_21, tau_a=tau_a)

        return np.concatenate([z_res_ary, z_ary_homo[0]]), np.concatenate([P_ary, P_ary_homo[0]]), ((dn_A_domega + dn_CMB_domega) / dn_CMB_domega), eps_max, omega_0
