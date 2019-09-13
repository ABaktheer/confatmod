'''
Created on 12.12.2016

@author: abaktheer
'''
from math import *

from matplotlib.figure import Figure
from scipy.interpolate import interp1d
from traits.api import \
    HasTraits, Property, Instance, cached_property, Str, Button, Enum, \
    Range, on_trait_change, Array, List, Float
from traitsui.api import \
    View, Item, Group, VGroup, HSplit, TreeEditor, TreeNode
from util.traits.editors import MPLFigureEditor
from view.ui import BMCSTreeNode, BMCSLeafNode

import matplotlib.gridspec as gridspec
from mats_bondslip_fatigue import MATSEvalFatigue
import numpy as np


class Material(BMCSLeafNode):

    node_name = Str('material parameters')

    E_b = Float(12900,
                label="E_b ",
                desc="Bond Stiffness",
                enter_set=True,
                auto_set=False)

    gamma = Float(55,
                  label="Gamma ",
                  desc="Kinematic hardening modulus",
                  enter_set=True,
                  auto_set=False)

    K = Float(11,
              label="K ",
              desc="Isotropic harening",
              enter_set=True,
              auto_set=False)

    S = Float(0.00051,
              label="S ",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    r = Float(0.51,
              label="r ",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    c = Float(2.8,
              label="c ",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    tau_pi_bar = Float(4.2,
                       label="Tau_pi_bar ",
                       desc="Reversibility limit",
                       enter_set=True,
                       auto_set=False)

    pressure = Float(0,
                     label="Pressure",
                     desc="Lateral pressure",
                     enter_set=True,
                     auto_set=False)

    a = Float(1.7,
              label="a",
              desc="Lateral pressure coefficient",
              enter_set=True,
              auto_set=False)

    view = View(VGroup(Group(Item('E_b'),
                             Item('tau_pi_bar'), show_border=True, label='Bond Stiffness and reversibility limit'),
                       Group(Item('gamma'),
                             Item('K'), show_border=True, label='Hardening parameters'),
                       Group(Item('S'),
                             Item('r'), Item('c'), show_border=True, label='Damage cumulation parameters'),
                       Group(Item('pressure'),
                             Item('a'), show_border=True, label='Lateral Pressure')))


class LoadingScenario(BMCSLeafNode):

    node_name = Str('Loading Scenario')
    number_of_cycles = Float(1.0)
    maximum_slip = Float(2)

    loading_type = Enum("Monotonic", "Cyclic")
    amplitude_type = Enum("Increased_Amplitude", "Constant_Amplitude")
    loading_range = Enum("Non_symmetric", "Symmetric")

    d_t = Float(0.005)
    t_max = Float(2.)
    k_max = Float(100)
    tolerance = Float(1e-4)

    d_array = Property(
        depends_on=' maximum_slip , number_of_cycles , loading_type , loading_range , amplitude_type , d_t , t_max')

    @cached_property
    def _get_d_array(self):
        number_of_increments = self.t_max / self.d_t
        if self.loading_type == "Monotonic":
            self.number_of_cycles = 1
        d_levels = np.linspace(0, self.maximum_slip, self.number_of_cycles * 2)
        d_levels[0] = 0

        if self.amplitude_type == "Increased_Amplitude" and self.loading_range == "Symmetric":
            d_levels.reshape(-1, 2)[:, 0] *= -1
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.amplitude_type == "Increased_Amplitude" and self.loading_range == "Non_symmetric":
            d_levels.reshape(-1, 2)[:, 0] *= 0
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.amplitude_type == "Constant_Amplitude" and self.loading_range == "Symmetric":
            d_levels.reshape(-1, 2)[:, 0] = -self.maximum_slip
            d_levels[0] = 0
            d_levels.reshape(-1, 2)[:, 1] = self.maximum_slip
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.amplitude_type == "Constant_Amplitude" and self.loading_range == "Non_symmetric":
            d_levels.reshape(-1, 2)[:, 0] *= 0
            d_levels.reshape(-1, 2)[:, 1] = self.maximum_slip
            s_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

    time_func = Property(depends_on='maximum_slip, t_max , d_array ')

    @cached_property
    def _get_time_func(self):
        t_arr = np.linspace(0, self.t_max, len(self.d_array))
        return interp1d(t_arr, self.d_array)

    figure = Instance(Figure)

    def _figure_default(self):
        figure = Figure()
        return figure

    update = Button()

    def _update_fired(self):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        x = np.arange(0, self.t_max, self.d_t)
        ax.plot(x, self.time_func(x))
        ax.set_xlabel('time')
        ax.set_ylabel('displacement')
        self.figure.canvas.draw()

    view = View(VGroup(Group(Item('loading_type')),
                       Group(Item('maximum_slip')),
                       Group(Item('number_of_cycles'),
                             Item('amplitude_type'),
                             Item('loading_range'), show_border=True, label='Cyclic load inputs'),
                       Group(Item('d_t'),
                             Item('t_max'),
                             Item('k_max'), show_border=True, label='Solver Settings')),
                Group(Item('update', label='Plot Loading senario')),
                Item('figure', editor=MPLFigureEditor(),
                     dock='horizontal', show_label=False))


class BondSlipModel(BMCSTreeNode):

    node_name = Str('Bond_slip_model')

    tree_node_list = List([])

    def _tree_node_list_default(self):
        # print 'NODE', self.material
        return [self.material, self.loading_scenario]

    material = Instance(Material)

    def _material_default(self):
        return Material()

    loading_scenario = Instance(LoadingScenario)

    def _loading_scenario_default(self):
        return LoadingScenario()

    mats_eval = Instance(MATSEvalFatigue)

    figure = Instance(Figure)

    def _figure_default(self):
        figure = Figure(facecolor='white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    def plot(self, figure, color='blue', linestyle='-',
             linewidth=1, label='<unnamed>'):
        pass
        # assign the material parameters
#         self.mats_eval.E_b = self.material.E_b
#         self.mats_eval.gamma = self.material.gamma
#         self.mats_eval.S = self.material.S
#         self.mats_eval.tau_pi_bar = self.material.tau_pi_bar
#         self.mats_eval.r = self.material.r
#         self.mats_eval.K = self.material.K
#         self.mats_eval.c = self.material.c
#         self.mats_eval.a = self.material.a
#         self.mats_eval.pressure = self.material.pressure
#
#         s_arr = self.loading_scenario._get_d_array()
#         t_arr = np.linspace(0, self.loading_scenario.t_max, len(s_arr))
#
#         tau_arr, w_arr, xs_pi_arr, xs_pi_cum = self.mats_eval.get_bond_slip(
#             s_arr)
#
#         ax1 = figure.add_subplot(221)
#         ax1.cla()
#         ax1.plot(s_arr, tau_arr, lw=linewidth, color=color,
#                  ls=linestyle, label=label)
#         ax1.set_title('Slip - Stress')
#         ax1.set_xlabel('Slip')
#         ax1.set_ylabel('Stress')
#         ax1.legend(loc=4)
#
#         ax2 = figure.add_subplot(222)
#         ax2.cla()
#         ax2.plot(s_arr, w_arr, lw=linewidth, color=color,
#                  ls=linestyle, label=label)
#         ax2.set_title('Slip - Damage')
#         ax2.set_xlabel('Slip')
#         ax2.set_ylabel('Damage')
#         ax2.set_ylim([0, 1])
#         ax2.legend(loc=4)
#
#         gs = gridspec.GridSpec(2, 2)
#         ax3 = figure.add_subplot(gs[-1, :])
#         ax3.plot(xs_pi_cum, w_arr, lw=linewidth, color=color,
#                  ls=linestyle, label=label)
#         ax3.set_title('Cumulative sliding - Damage')
#         ax3.set_ylim([0, 1])
#         ax3.set_xlabel('Cumulative sliding')
#         ax3.set_ylabel('Damage')
#         ax3.legend(loc=4)

    def plot_custom(self, ax1, ax2, ax3, ax4, ax5, ax6, color='blue', linestyle='-',
                    linewidth=1, alpha=1, label='<unnamed>'):

        self.mats_eval.E_b = self.material.E_b
        self.mats_eval.gamma = self.material.gamma
        self.mats_eval.S = self.material.S
        self.mats_eval.tau_pi_bar = self.material.tau_pi_bar
        self.mats_eval.r = self.material.r
        self.mats_eval.K = self.material.K
        self.mats_eval.c = self.material.c
        self.mats_eval.a = self.material.a
        self.mats_eval.pressure = self.material.pressure
        s_arr = self.loading_scenario._get_d_array()
        t_arr = np.linspace(0, self.loading_scenario.t_max, len(s_arr))

        tau_arr, w_arr, xs_pi_arr, xs_pi_cum, dissip, dissip_pi, dissip_w = self.mats_eval.get_bond_slip(
            s_arr)

        #======================================================================
        # Saving results
        #======================================================================
        n = len(s_arr)
        s_record = np.zeros(1)
        tau_record = np.zeros(1)
        w_record = np.zeros(1)
        xs_pi_record = np.zeros(1)
        xs_pi_cum_record = np.zeros(1)
        dissip_record = np.zeros(1)

        for i in range(0, n, 1):
            s_record = np.vstack((s_record, s_arr[i]))
            tau_record = np.vstack((tau_record, tau_arr[i]))
            #w_record = np.vstack((w_record, w_arr[i]))
            #xs_pi_record = np.vstack((xs_pi_record, xs_pi_arr[i]))
            #xs_pi_cum_record = np.vstack((xs_pi_cum_record, xs_pi_cum[i]))
            dissip_record = np.vstack((dissip_record, dissip[i]))

        np.savetxt(r'E:\Models_Implementation\python_results\bond_slip\slip.txt', np.transpose(
            s_record), delimiter=" ", fmt="%s")
        np.savetxt(r'E:\Models_Implementation\python_results\bond_slip\tau.txt', np.transpose(
            tau_record), delimiter=" ", fmt="%s")
#         np.savetxt(r'E:\Models_Implementation\python_results\bond_slip\w.txt', np.transpose(
#             w_record), delimiter=" ", fmt="%s")
#         np.savetxt(r'E:\Models_Implementation\python_results\bond_slip\xs_pi.txt', np.transpose(
#             xs_pi_record), delimiter=" ", fmt="%s")
#         np.savetxt(r'E:\Models_Implementation\python_results\bond_slip\xs_pi_cum.txt', np.transpose(
#             xs_pi_cum_record), delimiter=" ", fmt="%s")
        np.savetxt(r'E:\Models_Implementation\python_results\bond_slip\dissip.txt', np.transpose(
            dissip_record), delimiter=" ", fmt="%s")

        #================================================
        # plot stress-slip
        #================================================
        ax1.plot(s_arr, tau_arr, lw=linewidth, color=color,
                 ls=linestyle, alpha=alpha, label=label)
        #ax1.set_title('Slip - Stress')
        ax1.set_xlabel('Slip[mm]')
        ax1.set_ylabel('Stress[MPa]')
        ax1.legend(loc=4)

        #================================================
        # plot damage-slip
        #================================================
        ax2.plot(s_arr, w_arr, lw=linewidth, color=color,
                 ls=linestyle, alpha=alpha, label=label)
        #ax2.set_title('Slip - Damage')
        ax2.set_xlabel('Slip[mm]')
        ax2.set_ylabel('Damage')
        ax2.set_ylim([0, 1])
        ax2.legend(loc=4)

        #================================================
        # plot damage-cumulative sliding
        #================================================

        ax3.plot(xs_pi_cum, w_arr, lw=linewidth, color=color,
                 ls=linestyle, alpha=alpha, label=label)
        #ax3.set_title('Cumulative sliding - Damage')
        ax3.set_ylim([0, 1])
        ax3.set_xlabel('Cumulative sliding[mm]')
        ax3.set_ylabel('Damage')
        ax3.legend(loc=4)

        #================================================
        # plot slip-cumulative sliding
        #================================================
        ax4.plot(s_arr, xs_pi_cum, lw=linewidth, color=color,
                 ls=linestyle,  alpha=alpha, label=label)
        ax4.set_title('slip - Cumulative sliding')
        #ax4.set_ylim([0, 1])
        ax4.set_xlabel('slip[mm]')
        ax4.set_ylabel('Cumulative sliding[mm]')
        ax4.legend()

        #================================================
        # plot dissipation-damage
        #================================================
        ax5.plot(dissip, w_arr, lw=linewidth, color=color,
                 ls=linestyle, label=label)
        ax5.set_title('dissipation - damage')
        ax5.set_ylim([0, 1])
        ax5.set_xlabel('Dissipation')
        ax5.set_ylabel('Damage')
        ax5.legend()

        #================================================
        # plot dissipation-slip
        #================================================
        ax6.plot(s_arr, dissip, lw=linewidth, color=color,
                 ls=linestyle,  alpha=alpha, label=label)
        ax6.plot(s_arr, dissip_pi, lw=linewidth, color=color,
                 ls='--',  alpha=alpha, label=label)
        ax6.plot(s_arr, dissip_w, lw=linewidth, color=color,
                 ls=':',  alpha=alpha, label=label)
        #ax6.set_title('slip - dissipation')
        #ax6.set_ylim([0, 1])
        ax6.set_xlabel('Slip[mm]')
        ax6.set_ylabel('Dissipation[MPa]')
        # ax6.legend(loc=4)

#         #================================================
#         # plot time-stress
#         #================================================
#         ax5.plot(t_arr, tau_arr, lw=linewidth, color=color,
#                  ls=linestyle, alpha=alpha, label=label)
#         ax5.set_title('dissipation - damage')
#         #ax5.set_ylim([0, 1])
#         ax5.set_xlabel('time')
#         ax5.set_ylabel('stress[MPa]')
#         ax5.legend()


'''
if __name__ == '__main__':

    window = BondSlipModel(mats_eval=MATSEvalFatigue(), loading_scenario=LoadingScenario())
#     window.draw()
#
    window.configure_traits()
'''
