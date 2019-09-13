'''
Created on 12.01.2016
@author: ABaktheer
'''

#import math

#from ibvpy.api import BCDof
from matplotlib.figure import Figure
from scipy.interpolate import interp1d
from traits.api import \
    Property, Instance, cached_property, Str, Button, Enum, \
    Range, Array, List, Float, Int
from traitsui.api import \
    View, Item, Group, VGroup, HGroup
from util.traits.editors import MPLFigureEditor
from view.ui import BMCSTreeNode, BMCSLeafNode

import numpy as np

from .fets1d52ulrhfatigue import FETS1D52ULRHFatigue
from .mats_bondslip_fatigue import MATSEvalFatigue
from .tloop import TLoop
from .tstepper import TStepper


#from scipy import sparse
class Material(BMCSLeafNode):

    node_name = Str('material parameters')

    E_m = Float(30000,
                label="E_m ",
                desc="Matrix Stiffness",
                enter_set=True,
                auto_set=False)
    E_f = Float(230000,
                label="E_f ",
                desc="Fiber Stiffness",
                enter_set=True,
                auto_set=False)

    E_b = Float(12900,
                label="E_b ",
                desc="Bond Stiffness",
                enter_set=True,
                auto_set=False)

    gamma = Float(55.0,
                  label="Gamma ",
                  desc="Kinematic hardening modulus",
                  enter_set=True,
                  auto_set=False)

    K = Float(11.0,
              label="K ",
              desc="Isotropic harening",
              enter_set=True,
              auto_set=False)

    S = Float(0.00048,
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

    view = View(VGroup(Group(Item('E_m'),
                             Item('E_f'),
                             Item('E_b'),
                             Item('tau_pi_bar'), show_border=True, label='Bond Stiffness and reversibility limit'),
                       Group(Item('gamma'),
                             Item('K'), show_border=True, label='Hardening parameters'),
                       Group(Item('S'),
                             Item('r'), Item('c'), show_border=True, label='Damage cumulation parameters'),
                       Group(Item('pressure'),
                             Item('a'), show_border=True, label='Lateral Pressure')))


class Geometry(BMCSLeafNode):

    node_name = Str('Geometry')
    L_x = Range(1, 700, value=42)
    A_m = Float(15240.0, desc='matrix area [mm2]')
    A_f = Float(153.9, desc='reinforcement area [mm2]')
    P_b = Float(44.0, desc='perimeter of the bond interface [mm]')


class LoadingScenario(BMCSLeafNode):

    node_name = Str('Loading Scenario')
    number_of_cycles = Float(1.0)
    maximum_loading = Float(0.8)
    unloading_ratio = Range(0., 1., value=0.5)
    number_of_increments = Float(10)
    loading_type = Enum("Monotonic", "Cyclic")
    amplitude_type = Enum("Increased_Amplitude", "Constant_Amplitude")
    loading_range = Enum(
        "Non_symmetric", "Symmetric", "Loading blocks", "repeated loading blocks", "ten repeated loading blocks")

    cycles_n_1 = Float(10.0)
    cycles_n_2 = Float(10.0)
    cycles_n_3 = Float(10.0)
    cycles_n_4 = Float(10.0)
    cycles_n_5 = Float(10.0)

    maximum_loading_1 = Float(10.0)
    maximum_loading_2 = Float(10.0)
    maximum_loading_3 = Float(10.0)
    maximum_loading_4 = Float(10.0)
    maximum_loading_5 = Float(10.0)

    unloading_ratio_1 = Range(0., 1., value=0.1)
    unloading_ratio_2 = Range(0., 1., value=0.1)
    unloading_ratio_3 = Range(0., 1., value=0.1)
    unloading_ratio_4 = Range(0., 1., value=0.1)
    unloading_ratio_5 = Range(0., 1., value=0.1)

    max_monotonic = Float(10.0)
    number_of_levels = Int(10)
    number_of_cycles_each_level = Int(10)
    s_max_1 = Float(0.1)
    s_min = Float(0.05)
    s_max_2 = Float(0.7)

    number_of_repeted_blocks = Int(1)

    time = Range(0.00, 1.00, value=1.00)

    d_t = Float(0.005)
    t_max = Float(1.)
    k_max = Float(200)
    tolerance = Float(1e-4)

    d_array = Property(
        depends_on=' maximum_loading , number_of_cycles ,\
        loading_type , loading_range ,amplitude_type, unloading_ratio,\
        cycles_n_1, cycles_n_2,cycles_n_3,cycles_n_4,cycles_n_5,\
        maximum_loading_1,maximum_loading_2,maximum_loading_3,maximum_loading_4,maximum_loading_5,\
        unloading_ratio_1,unloading_ratio_2,unloading_ratio_3,unloading_ratio_4,unloading_ratio_5,\
        max_monotonic, number_of_levels,s_max_1 ,s_min ,s_max_2 , number_of_cycles_each_level, number_of_repeted_blocks')

    @cached_property
    def _get_d_array(self):

        if self.loading_type == "Monotonic":
            self.number_of_cycles = 1
            d_levels = np.linspace(
                0, self.maximum_loading, self.number_of_cycles * 2)
            d_levels[0] = 0
            d_levels.reshape(-1, 2)[:, 0] *= 0
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], self.number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.amplitude_type == "Increased_Amplitude" and self.loading_range == "Symmetric":
            d_levels = np.linspace(
                0, self.maximum_loading, self.number_of_cycles * 2)
            d_levels.reshape(-1, 2)[:, 0] *= -1
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], self.number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.amplitude_type == "Increased_Amplitude" and self.loading_range == "Non_symmetric":
            d_levels = np.linspace(
                0, self.maximum_loading, self.number_of_cycles * 2)
            d_levels.reshape(-1, 2)[:, 0] *= 0
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], self.number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.amplitude_type == "Constant_Amplitude" and self.loading_range == "Symmetric":
            d_levels = np.linspace(
                0, self.maximum_loading, self.number_of_cycles * 2)
            d_levels.reshape(-1, 2)[:, 0] = -self.maximum_loading
            d_levels[0] = 0
            d_levels.reshape(-1, 2)[:, 1] = self.maximum_loading
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], self.number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.amplitude_type == "Constant_Amplitude" and self.loading_range == "Non_symmetric":
            # d_1 = np.zeros(self.number_of_cycles*2 + 1)
            d_1 = np.zeros(1)
            d_2 = np.linspace(
                0, self.maximum_loading, self.number_of_cycles * 2)
            d_2.reshape(-1, 2)[:, 0] = self.maximum_loading
            d_2.reshape(-1, 2)[:, 1] = self.maximum_loading * \
                self.unloading_ratio
            d_history = d_2.flatten()
            d_arr = np.hstack((d_1, d_history))

            return d_arr

        if self.amplitude_type == "Constant_Amplitude" and self.loading_range == "Loading blocks":

            d_0 = np.zeros(1)

            d_1 = np.linspace(0, self.maximum_loading_1, self.cycles_n_1 * 2)
            d_1.reshape(-1, 2)[:, 0] = self.maximum_loading_1
            d_1.reshape(-1, 2)[:, 1] = self.maximum_loading_1 * \
                self.unloading_ratio_1
            d_history_1 = d_1.flatten()

            d_2 = np.linspace(0, self.maximum_loading_2, self.cycles_n_2 * 2)
            d_2.reshape(-1, 2)[:, 0] = self.maximum_loading_2
            d_2.reshape(-1, 2)[:, 1] = self.maximum_loading_2 * \
                self.unloading_ratio_2
            d_history_2 = d_2.flatten()

            d_arr = np.hstack((d_0, d_history_1, d_history_2))

            return d_arr

        if self.amplitude_type == "Constant_Amplitude" and self.loading_range == "repeated loading blocks":

            d_0 = np.zeros(1)

            d_1 = np.linspace(0, self.maximum_loading_1, self.cycles_n_1 * 2)
            d_1.reshape(-1, 2)[:, 0] = self.maximum_loading_1
            d_1.reshape(-1, 2)[:, 1] = self.maximum_loading_1 * \
                self.unloading_ratio_1
            d_history_1 = d_1.flatten()

            d_2 = np.linspace(0, self.maximum_loading_2, self.cycles_n_2 * 2)
            d_2.reshape(-1, 2)[:, 0] = self.maximum_loading_2
            d_2.reshape(-1, 2)[:, 1] = self.maximum_loading_2 * \
                self.unloading_ratio_2
            d_history_2 = d_2.flatten()

            d_3 = np.linspace(0, self.maximum_loading_3, self.cycles_n_3 * 2)
            d_3.reshape(-1, 2)[:, 0] = self.maximum_loading_3
            d_3.reshape(-1, 2)[:, 1] = self.maximum_loading_3 * \
                self.unloading_ratio_3
            d_history_3 = d_3.flatten()

            d_4 = np.linspace(0, self.maximum_loading_4, self.cycles_n_4 * 2)
            d_4.reshape(-1, 2)[:, 0] = self.maximum_loading_4
            d_4.reshape(-1, 2)[:, 1] = self.maximum_loading_4 * \
                self.unloading_ratio_4
            d_history_4 = d_4.flatten()

            d_5 = np.linspace(0, self.maximum_loading_5, self.cycles_n_5 * 2)
            d_5.reshape(-1, 2)[:, 0] = self.maximum_loading_5
            d_5.reshape(-1, 2)[:, 1] = self.maximum_loading_5 * \
                self.unloading_ratio_5
            d_history_5 = d_5.flatten()

            d_arr_1 = np.hstack(
                (d_0, d_history_1, d_history_2, d_history_3, d_history_4, d_history_5))

            d_arr = d_arr_1
            n = self.number_of_repeted_blocks
            for i in range(1, n):
                d_arr = np.hstack(
                    (d_arr, d_history_1, d_history_2, d_history_3, d_history_4, d_history_5))

            self.number_of_cycles = n * \
                (self.cycles_n_1 + self.cycles_n_2 +
                 self.cycles_n_3 + self.cycles_n_4 + self.cycles_n_5)
            return d_arr

        if self.amplitude_type == "Constant_Amplitude" and self.loading_range == "ten repeated loading blocks":
            d_0 = np.zeros(1)
            d_arr_1 = np.zeros(1)

            for i in range(1, self. number_of_levels + 1, 1):

                d_i = np.linspace(0, self.max_monotonic *
                                  (self.s_max_1 + (i - 1.0) * (self.s_max_2 -
                                                               self.s_max_1) / (self.number_of_levels - 1.0)),
                                  self.number_of_cycles_each_level * 2)
                d_i.reshape(-1, 2)[:, 0] = self.max_monotonic * (
                    self.s_max_1 + (i - 1.0) * (self.s_max_2 - self.s_max_1) / (self.number_of_levels - 1.0))
                d_i.reshape(-1, 2)[:, 1] = self.max_monotonic * self.s_min

                d_history_i = d_i.flatten()

                d_arr_1 = np.hstack((d_arr_1, d_history_i))

            d_arr = d_arr_1

            n = self.number_of_repeted_blocks
            for i in range(1, n):
                d_arr = np.hstack(
                    (d_arr, d_arr_1[1:]))

            self.number_of_cycles = n * \
                (self.number_of_levels * self.number_of_cycles_each_level)
            return d_arr

    time_func = Property(
        depends_on='maximum_loading, maximum_loading_1, maximum_loading_2,maximum_loading_3,maximum_loading_4,maximum_loading_5, t_max , d_array')

    @cached_property
    def _get_time_func(self):
        t_arr = np.linspace(0, self.t_max, len(self.d_array))
        print('time function'), interp1d(t_arr, self.d_array)
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
        ax.plot(x, self.time_func(x), 'k')
        ax.set_xlabel('time')
        ax.set_ylabel('displacement')
        self.figure.canvas.draw()

    view = View(VGroup(Group(Item('loading_type')),
                       Group(Item('maximum_loading')),
                       Group(Item('number_of_cycles'),
                             Item('amplitude_type'),
                             Item('loading_range'), Item('unloading_ratio'), show_border=True, label='Cyclic load inputs'),
                       Group(Item('d_t'),
                             Item('t_max'),
                             Item('k_max'), show_border=True, label='Solver Settings'),
                       Group(HGroup(Item('cycles_n_1'), Item('maximum_loading_1')), Item('unloading_ratio_1'),
                             HGroup(
                                 Item('cycles_n_2'), Item('maximum_loading_2')), Item('unloading_ratio_2'),
                             HGroup(
                                 Item('cycles_n_3'), Item('maximum_loading_3')), Item('unloading_ratio_3'),
                             HGroup(
                                 Item('cycles_n_4'), Item('maximum_loading_4')), Item('unloading_ratio_4'),
                             HGroup(
                                 Item('cycles_n_5'), Item('maximum_loading_5')), Item('unloading_ratio_5'),
                             Item('number_of_repeted_blocks'),
                             show_border=True, label='Loading Blocks')),
                HGroup(Item('max_monotonic'), Item('number_of_levels')), Item(
                    'number_of_cycles_each_level'), Item(
                    's_max_1'), Item('s_max_2'), Item('s_min'),
                Group(Item('update', label='Plot Loading scenario')),
                Item('figure', editor=MPLFigureEditor(),
                     dock='horizontal', show_label=False), Item('time', label='t/T_max'))


class PullOutSimulation(BMCSTreeNode):

    node_name = Str('pull out simulation')

    tree_node_list = List([])

    def _tree_node_list_default(self):

        return [self.material, self.geometry, self.loading_scenario]

    material = Instance(Material)

    def _material_default(self):
        return Material()

    loading_scenario = Instance(LoadingScenario)

    def _loading_scenario_default(self):
        return LoadingScenario()

    geometry = Instance(Geometry)

    def _geometry_default(self):
        return Geometry()

    mats_eval = Instance(MATSEvalFatigue)

    fets_eval = Instance(FETS1D52ULRHFatigue)

    time_stepper = Instance(TStepper)

    time_loop = Instance(TLoop)

    t_record = Array
    U_record = Array
    F_record = Array
    sf_record = Array
    eps_record = List
    sig_record = List
    D_record = List
    w_record = List

    figure = Instance(Figure)

    def _figure_default(self):
        figure = Figure()
        return figure

    def plot(self, figure, color='blue', linestyle='-',
             linewidth=1, label='<unnamed>'):
        pass

    def plot_custom(self, ax1, ax2, ax3, ax4, color='blue', linestyle='-',
                    linewidth=1, label='<unnamed>'):
        # assign the material parameters
        self.mats_eval.E_m = self.material.E_m
        self.mats_eval.E_f = self.material.E_f
        self.mats_eval.E_b = self.material.E_b
        self.mats_eval.gamma = self.material.gamma
        self.mats_eval.S = self.material.S
        self.mats_eval.tau_pi_bar = self.material.tau_pi_bar
        self.mats_eval.r = self.material.r
        self.mats_eval.K = self.material.K
        self.mats_eval.c = self.material.c
        self.mats_eval.a = self.material.a
        self.mats_eval.pressure = self.material.pressure
        # assign the geometry parameters
        self.fets_eval.A_m = self.geometry.A_m
        self.fets_eval.P_b = self.geometry.P_b
        self.fets_eval.A_f = self.geometry.A_f
        self.time_stepper.L_x = self.geometry.L_x
        # assign the parameters for solver and loading_scenario
        self.time_loop.t_max = self.loading_scenario.t_max
        self.time_loop.d_t = self.loading_scenario.d_t
        self.time_loop.k_max = self.loading_scenario.k_max
        self.time_loop.tolerance = self.loading_scenario.tolerance
        # assign the bc
        self.time_stepper.bc_list[1].value = 1
        self.time_stepper.bc_list[
            1].time_function = self.loading_scenario.time_func

        self.U_record, self.F_record, self.sf_record, self.t_record, self.eps_record, \
            self.sig_record, self.w_record, self.D_record = self.time_loop.eval()
        n_dof = 2 * self.time_stepper.domain.n_active_elems + 1

        #=======================================================
        # plotting the pull-out curve
        #=======================================================
        if self.loading_scenario.loading_type == "Monotonic" or "Cyclic":

            # Unloaded end
            ax1.plot(self.U_record[:, 1], self.F_record[:, n_dof] / 1000, lw=linewidth, color=color,
                     ls=linestyle, label=label)

            # loaded end
            ax1.plot(self.U_record[:, n_dof], self.F_record[:, n_dof] / 1000.0, lw=linewidth, color=color,
                     ls=linestyle, label=label)
            print(np.nanmax(self.F_record[:, n_dof]))
            # ax1.plot(self.U_record, self.F_record / 1000, lw=linewidth, color=color,
            #         ls=linestyle, label=label)

#             # Eligehausen test
#             slip, force = np.loadtxt(
#                 'E:\\publishing\Euro-C-2018\Pull_out_Test\Results\Rehm_pullout_fatigue\Monotonic-2.txt')  # Experimental results
#             ax1.plot(slip, force, '--k', label='Experiment')

#             # FRP (Carloni) test
#             slip, force = np.loadtxt(
#                 r'E:\\Models_Implementation\FRP_pullout_fatigue\results\Exp\monotonic.txt')
#             ax1.plot(slip, force, '--k', label='Experiment')

            #=======================
            # IMB - Beam End test
            #=======================

            # LS1 # ICEM
            u_1, f_1 = np.loadtxt(
                r'E:\Publishing\ICEM_2018\results\exp\monotonic_16_01.txt')
            u_2, f_2 = np.loadtxt(
                r'E:\Publishing\ICEM_2018\results\exp\monotonic_16_02.txt')
            u_3, f_3 = np.loadtxt(
                r'E:\Publishing\ICEM_2018\results\exp\monotonic_16_03.txt')
            ax1.plot(f_1, u_1, 'k', label='1')
            ax1.plot(f_2, u_2, 'k', label='2')
            ax1.plot(f_3, u_3, 'k', label='3')

#             # LS1 # FIB
#             u_1, f_1 = np.loadtxt(
#                 r'E:\Publishing\FIB_2018_(2)\results\exp\monotonic_16_active.txt')
#             u_2, f_2 = np.loadtxt(
#                 r'E:\Publishing\FIB_2018_(2)\results\exp\monotonic_16_passiv.txt')
#             ax1.plot(f_1, u_1, 'k', label='loaded')
#             ax1.plot(f_2, u_2, '--k', label='unloaded')


#             # LS2
#             u_1, f_1 = np.loadtxt(
#                 r'E:\Publishing\FIB_2018_(2)\results\LS2_F_w.txt')
#             ax1.plot(f_1, u_1, 'k', label='loaded')

            # for saving results
            n = len(self.U_record[:, n_dof])
            U_record_1 = np.zeros(1)
            U_record_2 = np.zeros(1)
            F_record = np.zeros(1)

            for i in range(0, n, 1):
                U_record_1 = np.vstack((U_record_1, self.U_record[i, n_dof]))
                U_record_2 = np.vstack((U_record_2, self.U_record[i, 1]))
                F_record = np.vstack(
                    (F_record, self.F_record[i, n_dof] / 1000.0))

            # saving results
            np.savetxt(r'e:\\Models_implementation\python_results\U_loaded_slip2.txt', np.transpose(
                U_record_1), delimiter=" ", fmt="%s")
            np.savetxt(r'e:\\Models_implementation\python_results\U_unloaded_slip2.txt', np.transpose(
                U_record_2), delimiter=" ", fmt="%s")
            np.savetxt(r'e:\\Models_implementation\python_results\F2.txt', np.transpose(
                F_record), delimiter=" ", fmt="%s")

            ax1.set_title('Pull-out curve')
            ax1.set_xlabel('Slip(mm)')
            ax1.set_ylabel('Force (KN)')
            #ax1.set_xlim(0, 1.5)
            # ax1.legend(loc=4)


#

        #================================================================
        # plotting the max slip for each cycle (S(n) curve) - logarithmic
        #================================================================
        n = (len(self.loading_scenario.d_array) - 1) / 2
        print('n', n)

        if self.loading_scenario.loading_type == "Cyclic":
            u_max_1 = np.zeros(1)
            u_max_2 = np.zeros(1)
            f_max_2 = np.zeros(1)
            f_min_2 = np.zeros(1)
            u_min_1 = np.zeros(1)
            u_min_2 = np.zeros(1)
            # E_ed = np.zeros(n)
            t = np.zeros(1)

            for i in range(0, n, 1):
                idx = (2 * i + 1) * (self.loading_scenario.t_max) / \
                    (2 * n * self.loading_scenario.d_t)
                idx = int(idx)
                print('idx', idx)
#                 idx_2 = ( i ) * (self.loading_scenario.t_max) / \
#                     (2 * n * self.loading_scenario.d_t)
                if idx >= len(self.t_record):
                    break
                else:
                    # max slip of the loaded end
                    print('u_max_1', u_max_1)
                    print('n_dof', n_dof)
                    print('self.U_record', self.U_record)
                    u_max_1 = np.vstack((u_max_1, self.U_record[idx, n_dof]))

                    # max slip of the unloaded end
                    u_max_2 = np.vstack((u_max_2, self.U_record[idx, 1]))
                    f_max_2 = np.vstack((f_max_2, self.F_record[idx, n_dof]))
                    #u_min_1 = np.vstack((u_min_1, self.U_record[idx_2, n_dof]))

                    t = np.vstack((t, self.t_record[idx]))

            for i in range(1, n, 1):

                idx_2 = (2 * i) * (self.loading_scenario.t_max) / \
                    (2 * n * self.loading_scenario.d_t)
                idx_2 = int(idx_2)
                if idx_2 >= len(self.t_record):
                    break
                else:
                    u_min_1 = np.vstack((u_min_1, self.U_record[idx_2, n_dof]))
                    u_min_2 = np.vstack((u_min_2, self.U_record[idx_2, 1]))
                    f_min_2 = np.vstack((f_min_2, self.F_record[idx_2, n_dof]))

        if self.loading_scenario.loading_type == "Cyclic":

            '''not logarithmic'''
            ax3.plot(t[1:-1] * (self.loading_scenario.number_of_cycles / self.loading_scenario.t_max), u_max_1[1:-1],
                     lw=linewidth, color=color, ls=linestyle, label=label)
            ax3.plot(t[1:-1] * (self.loading_scenario.number_of_cycles / self.loading_scenario.t_max), u_max_2[1:-1],
                     lw=linewidth, color=color, ls=linestyle, label=label)

            '''half logarithmic'''
#             ax3.plot(np.log10(t[1:-1] * (self.loading_scenario.number_of_cycles / self.loading_scenario.t_max)), u_max_2[1:-1],
#                      lw=linewidth, color=color, ls=linestyle, label=label)

            # saving results
            np.savetxt(r'e:\\Models_implementation\python_results\n2.txt', np.transpose(t[
                       0:-1] * (self.loading_scenario.number_of_cycles / self.loading_scenario.t_max)), delimiter=" ", fmt="%s")
            np.savetxt(r'e:\\Models_implementation\python_results\s_unloaded2.txt', np.transpose(
                u_max_2[0:-1]), delimiter=" ", fmt="%s")
            np.savetxt(r'e:\\Models_implementation\python_results\s_loaded2.txt', np.transpose(
                u_max_1[0:-1]), delimiter=" ", fmt="%s")
            np.savetxt(r'e:\\Models_implementation\python_results\s_min_loaded.txt', np.transpose(
                u_min_1[0:-1]), delimiter=" ", fmt="%s")
            np.savetxt(r'e:\\Models_implementation\python_results\s_min_unloaded.txt', np.transpose(
                u_min_2[0:-1]), delimiter=" ", fmt="%s")
            np.savetxt(r'e:\\Models_implementation\python_results\F_max_unloaded.txt', np.transpose(
                f_max_2[0:-1]), delimiter=" ", fmt="%s")
            np.savetxt(r'e:\\Models_implementation\python_results\F_min_unloaded.txt', np.transpose(
                f_min_2[0:-1]), delimiter=" ", fmt="%s")


#             ax3.plot(t[1:-1] * (self.loading_scenario.number_of_cycles / self.loading_scenario.t_max), u_max_2[1:-1],
#                      lw=linewidth, color=color, ls=linestyle, label=label)

            ax3.set_ylim(0.0, 2.0)

            ax3.set_title('S(n)')
            ax3.set_xlabel('Log(N)')
            ax3.set_ylabel('Max Slip(mm)')

# #         #================================================================
# #         # plotting the damage distribution along the bonded length
# #         #================================================================
#         X = np.linspace(0, self.time_stepper.L_x, self.time_stepper.n_e_x + 1)
#         X_ip = np.repeat(X, 2)[1:-1]
#         idx_arr = np.zeros(n)
#         for i in range(0, n, 1):
#             idx = (2 * i + 1) * (self.loading_scenario.t_max) / \
#                 (2 * n * self.loading_scenario.d_t)
#
#             if idx >= len(self.t_record):
#                 break
#             idx_arr[i] = idx
#
#         print 'idx_arr', idx_arr
#
# #         for i in [1, 5, 10, -1]:
# #
# #             ax2.plot(X_ip, self.w_record[int(idx_arr[i])].flatten(), lw=linewidth, color=color,
# #                      ls=linestyle, label=label)
#
#         np.savetxt(r'e:\\Models_implementation\python_results\interface.txt',
#                    X_ip, delimiter=" ", fmt="%s")
# #             np.savetxt(r'e:\\Models_implementation\python_results\w_cycle_1.txt',
# # self.w_record[int(idx_arr[1])].flatten().T, delimiter=" ", fmt="%s")
#
#         n = int(self.loading_scenario.t_max / self.loading_scenario.d_t)
#         w_arr = np.zeros(1)
#         idx = idx_arr.astype(int)
#         for i in idx:
#             if i <= len(self.t_record) - 1:
#                 w_arr = np.hstack(
#                     (w_arr, self.w_record[i].flatten().T))
#             else:
#                 break
#
#         np.savetxt(r'E:\Models_Implementation\python_results\damage.txt',
#                    w_arr, delimiter=" ", fmt="%s")

#         #================================================================
#         # plotting the strain distribution along the bonded length
#         #================================================================
#         X = np.linspace(0, self.time_stepper.L_x, self.time_stepper.n_e_x + 1)
#         X_ip = np.repeat(X, 2)[1:-1]
#         idx_arr = np.zeros(n)
# #         for i in range(0, n, 1):
# #             idx = (2 * i + 1) * (self.loading_scenario.t_max) / \
# #                 (2 * n * self.loading_scenario.d_t)
# #             # print 'idx', idx
# #             if idx >= len(self.t_record):
# #                 break
# #             # else:
# #             idx_arr[i] = idx
# #
# #         print 'idx_arr', idx_arr
#
#         # for i in range(0, len(self.loading_scenario.d_array), 1):
#
# #         ax2.plot(X_ip, self.eps_record[i][:, :, 0].flatten(), lw=linewidth, color=color,
# #                      ls=linestyle, label=label)
#
# #         ax2.plot(X_ip, self.eps_record[i - 1][:, :, 2].flatten(), lw=linewidth, color=color,
# #                  ls=linestyle, label=label)
#
#         #self.sig_record[-1][:, :, 0].flatten()
#
#         np.savetxt(r'e:\\Models_implementation\python_results\interface.txt',
#                    X_ip, delimiter=" ", fmt="%s")
#
#         n = int(self.loading_scenario.t_max / self.loading_scenario.d_t)
#
#         # print 'n', n
#         eps_arr = np.zeros(1)
#
#         for i in range(0, n, 1):
#             if i <= len(self.t_record) - 1:
#                 eps_arr = np.hstack(
#                     (eps_arr, self.eps_record[i][:, :, 2].flatten().T))
#             else:
#                 break
#
#         np.savetxt(r'E:\Models_Implementation\python_results\strain\fiber_strain.txt',
#                    eps_arr, delimiter=" ", fmt="%s")

        #======================================================================
        # plotting the damage evolution at different points along the bonded length
        #======================================================================
#         if self.loading_scenario.loading_type == "Monotonic":
#             n_e = self.time_stepper.domain.n_active_elems
#             # print 't', self.t_record.shape
#             # print 'w', self.w_record[:][:]
#             w_first = np.concatenate(self.w_record[:][:], axis=0)[0::n_e][:, 0]
#             w_last = np.concatenate(
#                 self.w_record[:][:], axis=0)[n_e - 1::n_e][:, -1]
#
#             # print 'w_first', w_first
#             # for saving results
#             m = len(w_first)
#             w_first_1 = np.zeros(1)
#             w_last_1 = np.zeros(1)
#             time = np.zeros(1)
#             for i in range(0, m, 1):
#                 w_first_1 = np.vstack((w_first_1, w_first[i]))
#                 w_last_1 = np.vstack((w_last_1, w_last[i]))
#                 time = np.vstack(
#                     (time, self.t_record[i] / self.t_record[-1]))
#
#             # saving results
#             np.savetxt(r'e:\\Models_implementation\python_results\M_time.txt', np.transpose(
#                 time), delimiter=" ", fmt="%s")
#             np.savetxt(r'e:\\Models_implementation\python_results\M_w_Unloaded.txt', np.transpose(
#                 w_first_1), delimiter=" ", fmt="%s")
#             np.savetxt(r'e:\\Models_implementation\python_results\M_w_Loaded.txt', np.transpose(
#                 w_last_1), delimiter=" ", fmt="%s")
#
#             ax2.plot(self.t_record / self.t_record[-1], w_first, lw=linewidth, color=color,
#                      ls='--', label='unloaded end')
#             ax2.plot(self.t_record / self.t_record[-1], w_last, lw=linewidth, color=color,
#                      ls='-', label='loaded end')
#
#             ax2.set_title('damage evolution')
#             ax2.set_xlabel('time')
#             ax2.set_ylabel('damage')
#             ax2.set_ylim(0, 1.0)
#             ax2.legend(loc=4)
#
#         if self.loading_scenario.loading_type == "Cyclic":
#             n_e = self.time_stepper.domain.n_active_elems
#             n = (len(self.loading_scenario.d_array) - 1) / 2
#
#             w_first = np.concatenate(self.w_record[:][:], axis=0)[0::n_e][:, 0]
#             w_last = np.concatenate(
#                 self.w_record[:][:], axis=0)[n_e - 1::n_e][:, -1]
#
#             w_first_cycle = np.zeros(1)
#             w_last_cycle = np.zeros(1)
#             t_cycle = np.zeros(1)
#
#             for i in range(0, n, 1):
#                 idx = (2 * i + 1) * (self.loading_scenario.t_max) / \
#                     (2 * n * self.loading_scenario.d_t)
#                 if idx >= len(self.t_record):
#                     break
#                 else:
#                     w_first_cycle = np.vstack((w_first_cycle, w_first[idx]))
#                     w_last_cycle = np.vstack((w_last_cycle, w_last[idx]))
#                     t_cycle = np.vstack((t_cycle, self.t_record[idx]))
#
#             w_first_cycle[-2] = 1.0
#             w_last_cycle[-2] = 1.0
#
#             # saving results
#             np.savetxt(r'e:\\Models_implementation\python_results\C_time.txt', np.transpose(
#                 t_cycle[0:-1] / t_cycle[-1]), delimiter=" ", fmt="%s")
#             np.savetxt(r'e:\\Models_implementation\python_results\C_w_Unloaded.txt', np.transpose(
#                 w_first_cycle[0:-1]), delimiter=" ", fmt="%s")
#             np.savetxt(r'e:\\Models_implementation\python_results\C_w_Loaded.txt', np.transpose(
#                 w_last_cycle[0:-1]), delimiter=" ", fmt="%s")
#
#             ax2.plot(t_cycle[0:-1] / t_cycle[-1],  w_first_cycle[0:-1], lw=linewidth, color=color,
#                      ls='--', label='unloaded end')
#             ax2.plot(t_cycle[0:-1] / t_cycle[-1], w_last_cycle[0:-1], lw=linewidth, color=color,
#                      ls='-', label='loaded end')
#
#             ax2.set_title('damage evolution')
#             ax2.set_xlabel('time')
#             ax2.set_ylabel('damage')
#             ax2.set_ylim(0, 1.0)
#             ax2.legend(loc=4)


#         #======================================================================
#         # plotting the max slip for each cycle (S(n) curve) - logarithmic up to delta s= 0.1 mm
#         #======================================================================
#         n = (len(self.loading_scenario.d_array) - 1) / 2
#         u_max_1 = np.zeros(1)
#         u_max_2 = np.zeros(1)
#         u_min = np.zeros(1)
#         # E_ed = np.zeros(n)
#         t = np.zeros(1)
#
#         for i in range(0, n, 1):
#             idx = (2 * i + 1) * (self.loading_scenario.t_max) / \
#                 (2 * n * self.loading_scenario.d_t)
#             if idx >= len(self.t_record):
#                 break
#             else:
#                 # max slip of the loaded end
#                 #u_max_1 = np.vstack((u_max_1, self.U_record[idx, n_dof]))
#                 # max slip of the unloaded end
#                 u_max_2 = np.vstack((u_max_2, self.U_record[idx, 1]))
#                 t = np.vstack((t, self.t_record[idx]))
#
# #                 if u_max_2[-1] - u_max_2[1] >= 0.1:
# #                     break
#
#         if self.loading_scenario.loading_type == "Cyclic":
#             #             ax4.plot(t[1:-1] * (self.loading_scenario.number_of_cycles / self.loading_scenario.t_max), u_max_2[1:-1],
#             # lw=linewidth, color=color, ls=linestyle, label=label)
#
#             ax4.plot(t[1:-1] / t[-1], u_max_2[1:-1],
#                      lw=linewidth, color=color, ls=linestyle, label=label)
#
#             # ax4.set_xscale("log")
#             ax4.set_ylim(0, 2.5)
#             #ax4.set_xlim(0, 6.0)
#             # ax4.set_xlim(left=1.0)
#
#             ax4.set_title('Creep-fatigue curve')
#             ax4.set_xlabel('Number of cycles(n)')
#             # ax4.set_xlabel('N/Nf')
#             ax4.set_ylabel('Max Slip(mm)')
#             ax4.legend(loc=2)

        #======================================================
        # plotting the max slip for each cycle (S(n) curve) double logarithmic
        #======================================================
#             ax4.loglog(t[1:-1] * (self.loading_scenario.number_of_cycles / self.loading_scenario.t_max), u_max_2[1:-1],
#                        lw=linewidth, color=color, ls=linestyle, label=label)
#
#             #ax4.plot(N_1, s_1, '--k', label='Experiment, S=0.85')
#             #ax4.plot(N_2, s_2, '--k', label='Experiment, S=0.77')
#             ax4.loglog(N_1, s_1, '--k')  # , label='Experiment, S=0.85')
#             ax4.loglog(N_2, s_2, '--k')  # , label='Experiment, S=0.77')
# #             ax4.loglog(N_3, s_3, '--k')  # , label='Experiment, S=0.65')
# #             ax4.loglog(N_4, s_4, '--k')  # , label='Experiment, S=0.50')
# #             ax4.loglog(N_5, s_5, '--k')  # , label='Experiment, S=0.40')
#             ax4.set_xscale("log")
#             ax4.set_yscale("log")
#             ax4.set_xlim(1e0, 1e6)
#             ax4.set_ylim(1e-1, 1.5)
#             ax4.set_title('S(n)')
#             ax4.set_xlabel('Number of cycles(n)')
#             ax4.set_ylabel('Max Slip(mm)')
#             # ax4.legend(loc=4)
#

        #==============================================
        # plotting the stiffness vs. number of cycles
        #==============================================

#         if   self.loading_scenario.loading_type == "Cyclic":
#             ax3.plot(t[1:-1] *(self.loading_scenario.number_of_cycles / self.loading_scenario.t_max),
#                       (self.loading_scenario.maximum_loading - self.loading_scenario.maximum_loading *
#                        self.loading_scenario.unloading_ratio) /
#                                 (u_max[1:-1]),
#                                 lw=linewidth, color=color,
#                                 ls=linestyle, label=label)
#             ax3.set_xlim(0, 1)
#             ax3.set_title('Stiffness vs.  number of cycles')
#             ax3.set_xlabel('N'); ax3.set_ylabel('Stiffness')
#             ax3.legend(loc=4)

        #==============================================
        # Extra plots - plot only specific cycles
        #==============================================

#         if self.loading_scenario.loading_type == "Cyclic":
#
#             r = 3
#             m = (2 * r + 1) * ((self.loading_scenario.t_max) /
#                                (2 * n * self.loading_scenario.d_t)) + r
#             U_sp = np.zeros(m)
#             F_sp = np.zeros(m)
#
#             idx_0 = (self.loading_scenario.t_max) / \
#                 (2 * n * self.loading_scenario.d_t) + 1
#             U_sp[:idx_0] = self.U_record[:idx_0, n_dof]
#             F_sp[:idx_0] = self.F_record[:idx_0, n_dof]
#
#             idx_1 = (2 * 1 + 1) * (self.loading_scenario.t_max) / \
#                 (2 * n * self.loading_scenario.d_t) + 1
#             U_sp[idx_0:idx_1] = self.U_record[idx_0:idx_1, n_dof]
#             F_sp[idx_0:idx_1] = self.F_record[idx_0:idx_1, n_dof]
#
#             idx_2 = (2 * 400 + 1) * (self.loading_scenario.t_max) / \
#                 (2 * n * self.loading_scenario.d_t)
#             idx_22 = (2 * 401 + 1) * (self.loading_scenario.t_max) / \
#                 (2 * n * self.loading_scenario.d_t) + 1
#             U_sp[idx_1:idx_1 + 2 * idx_0 -
#                  1] = self.U_record[idx_2:idx_22, n_dof]
#             F_sp[idx_1:idx_1 + 2 * idx_0 -
#                  1] = self.F_record[idx_2:idx_22, n_dof]
#
#             idx_3 = (2 * 828 + 1) * (self.loading_scenario.t_max) / \
#                 (2 * n * self.loading_scenario.d_t)
#             idx_33 = (2 * 829 + 1) * (self.loading_scenario.t_max) / \
#                 (2 * n * self.loading_scenario.d_t) + 1
#
#             U_sp[idx_1 + 2 * idx_0 - 1: idx_1 + 4 * idx_0 -
#                  2] = self.U_record[idx_3: idx_33, n_dof]
#             F_sp[idx_1 + 2 * idx_0 - 1: idx_1 + 4 * idx_0 -
#                  2] = self.F_record[idx_3: idx_33, n_dof]
#
#             ax1.plot(U_sp, F_sp / 1000, lw=linewidth,
#                      color=color, ls=linestyle, label=label)
#             ax1.set_title('Pull-out curve')
#             ax1.set_xlabel('Slip(mm)')
#             ax1.set_ylabel('Force (KN)')
#             ax1.legend(loc=4)
