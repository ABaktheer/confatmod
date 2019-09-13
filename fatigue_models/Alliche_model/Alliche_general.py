'''
Created on 14.06.2017

@author: abaktheer
'''

'''
Implementation of the fatigue model for plain concrete [A.Alliche, 2004] under uniaxial compressive loading
(stress driven algorithm)
'''


from matplotlib.figure import Figure
from scipy.interpolate import interp1d
from traits.api import Property, Instance,\
    cached_property, Str, Button, Enum, Float, List, Range, Int
from traitsui.api import View, Item, Group, VGroup, HGroup, TreeNode
from util.traits.editors import MPLFigureEditor
from view.ui import BMCSLeafNode, BMCSTreeNode
from view.window import BMCSWindow
import numpy as np


class Material(BMCSLeafNode):
    #--------------------------
    # model parameters
    #--------------------------

    node_name = Str('Material')

    lamda = Float(13972.2,
                  label="lamda",
                  desc="first lame constant",
                  enter_set=True,
                  auto_set=False)

    mu = Float(20958.3,
               label="mu",
               desc="second lame constant",
               enter_set=True,
               auto_set=False)

    alpha = Float(2237.5,
                  label="alpha",
                  desc="---",
                  enter_set=True,
                  auto_set=False)

    beta = Float(-2116.5,
                 label="beta",
                 desc="---",
                 enter_set=True,
                 auto_set=False)

    g = Float(-10.0,
              label="g",
              desc="---",
              enter_set=True,
              auto_set=False)

    C0 = Float(0.0,
               label="C0",
               desc="---",
               enter_set=True,
               auto_set=False)

    C1 = Float(0.0019,
               label="C1",
               desc="---",
               enter_set=True,
               auto_set=False)

    K = Float(0.003345,
              label="K",
              desc="---",
              enter_set=True,
              auto_set=False)

    n = Float(10,
              label="n",
              desc="---",
              enter_set=True,
              auto_set=False)


class MATSEvalAllicheFatigue(BMCSLeafNode):
    #--------------------------
    # model parameters
    #--------------------------
    lamda = Float(13972.2,
                  label="lamda",
                  desc="first lame constant",
                  enter_set=True,
                  auto_set=False)

    mu = Float(20958.3,
               label="mu",
               desc="second lame constant",
               enter_set=True,
               auto_set=False)

    alpha = Float(2237.5,
                  label="alpha",
                  desc="---",
                  enter_set=True,
                  auto_set=False)

    beta = Float(-2116.5,
                 label="beta",
                 desc="---",
                 enter_set=True,
                 auto_set=False)

    g = Float(-9.9,
              label="g",
              desc="---",
              enter_set=True,
              auto_set=False)

    C0 = Float(0.0,
               label="C0",
               desc="---",
               enter_set=True,
               auto_set=False)

    C1 = Float(0.002033,
               label="C1",
               desc="---",
               enter_set=True,
               auto_set=False)

    K = Float(0.003345,
              label="K",
              desc="---",
              enter_set=True,
              auto_set=False)

    n = Float(10,
              label="n",
              desc="---",
              enter_set=True,
              auto_set=False)

    def get_stress_strain(self, sigma_1_arr):

        #----------------------------------------------------------------------
        # arrays to store the values
        #----------------------------------------------------------------------
        # normal strain
        eps_1_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
        # lateral strain
        eps_2_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
        # damage factor
        w_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
        f_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)

    #-----------------------------------------------------------------------
    # material parameters
    #-----------------------------------------------------------------------
        # lame constants [MPa]
        lamda = self.lamda
        mu = self.mu
        # fatigue model material parameter
        alpha = self.alpha
        beta = self.beta
        g = self.g
        C0 = self.C0
        C1 = self.C1
        K = self.K
        n = self.n

    #-----------------------------------------------------------------------
    # state variables
    #-----------------------------------------------------------------------
        #sigma_1_arr[0] = 0
        eps_1_i = 0.0
        eps_2_i = 0.0
        w_i = 0.0

        for i in range(1, len(sigma_1_arr)):
            print(i)
            sigma_1_i = sigma_1_arr[i]

            eps_2_i = -1.0 * ((lamda + alpha * w_i) *
                              sigma_1_i + g * w_i * (lamda + 2.0 * mu)) / \
                ((lamda + 2.0 * mu) * (2.0 * (lamda + mu) + 4.0 *
                                       w_i * (alpha + beta)) - 2.0 * (lamda + alpha * w_i) ** 2)

            eps_1_i = sigma_1_i / \
                (lamda + 2.0 * mu) - 2.0 * eps_2_i *  \
                (lamda + alpha * w_i) / (lamda + 2.0 * mu)

            f_i = abs(g) * eps_2_i - (C0 + 2 * C1 * w_i)

            kappa_i = (lamda + 2.0 * mu) * (2.0 * (lamda + mu) +
                                            4.0 * w_i * (alpha + beta) -
                                            alpha * (g / (2.0 * C1)) *
                                            (2.0 * eps_2_i + eps_1_i) -
                                            (g**2.0 / (2.0 * C1))) - 2.0 * (lamda + alpha * w_i)**2

            d_sigma_1 = sigma_1_arr[i] - sigma_1_arr[i - 1]
            m = -1.0 * ((lamda + alpha * w_i) / kappa_i) * d_sigma_1

            # loading stage (evolve of the fatigue damage based on (Marigo.85)
            # model)
            if m > 0:
                d_w = m * abs(g) / (2.0 * C1) * (f_i / K)**n
            else:  # unloading stage (no fatigue damage)
                d_w = 0

            w_i = w_i + d_w

            if w_i > 5.0:
                print(' ----------> No Convergence any more')
                break

            if abs(eps_1_i) > 0.0035:
                print(' ----------> No Convergence any more')
                break

            eps_1_arr[i] = eps_1_i
            eps_2_arr[i] = eps_2_i
            w_arr[i] = w_i
            f_arr[i] = f_i

        return sigma_1_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, i


class LoadingScenario(BMCSLeafNode):

    node_name = Str('Loading Scenario')
    number_of_cycles = Float(1.0)
    maximum_loading = Float(0.8)
    unloading_ratio = Range(0., 1., value=0.5)

    loading_type = Enum("Monotonic", "Cyclic")
    amplitude_type = Enum("Constant_Amplitude", "Variable_Amplitude")
    loading_levels = Enum(
        "Increased_blocks(LS2)", "Flexible_loading_levels(LS3,LS4)")

    cycles_n_1 = Float(10.0)
    cycles_n_2 = Float(10.0)
    cycles_n_3 = Float(10.0)
    cycles_n_4 = Float(10.0)
    cycles_n_5 = Float(10.0)

    S_max_1 = Float(1.0)
    S_max_2 = Float(1.0)
    S_max_3 = Float(1.0)
    S_max_4 = Float(1.0)
    S_max_5 = Float(1.0)

    S_min_1 = Range(0., 1., value=0.1)
    S_min_2 = Range(0., 1., value=0.1)
    S_min_3 = Range(0., 1., value=0.1)
    S_min_4 = Range(0., 1., value=0.1)
    S_min_5 = Range(0., 1., value=0.1)

    max_monotonic = Float(10.0)
    number_of_levels = Int(10)
    number_of_cycles_each_level = Int(10)
    s_max_first = Float(0.1)
    s_min_all = Float(0.05)
    s_max_last = Float(0.7)

    number_of_repeted_blocks = Int(2)

    time = Range(0.00, 1.00, value=1.00)

    d_t = Float(0.005)
    t_max = Float(1.)

    d_array = Property(
        depends_on=' d_t, t_max, maximum_loading , number_of_cycles ,\
        loading_type ,amplitude_type, loading_levels , unloading_ratio,\
        cycles_n_1, cycles_n_2,cycles_n_3,cycles_n_4,cycles_n_5,\
        S_max_1,S_max_2,S_max_3,S_max_4,S_max_5,\
        S_min_1,S_min_2,S_min_3,S_min_4,S_min_5,\
        max_monotonic, number_of_levels,s_max_first ,s_min_all ,s_max_last , number_of_cycles_each_level, number_of_repeted_blocks')

    @cached_property
    def _get_d_array(self):
        number_of_increments = self.t_max / self.d_t
        if self.loading_type == "Monotonic":
            self.number_of_cycles = 1
            d_levels = np.linspace(
                0, self.maximum_loading, self.number_of_cycles * 2)
            d_levels[0] = 0
            d_levels.reshape(-1, 2)[:, 0] *= 0
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.loading_type == "Cyclic" and self.amplitude_type == "Constant_Amplitude":
            d_1 = np.zeros(1)
            d_2 = np.linspace(
                0, self.maximum_loading, self.number_of_cycles * 2)
            d_2.reshape(-1, 2)[:, 0] = self.maximum_loading
            d_2.reshape(-1, 2)[:, 1] = self.maximum_loading * \
                self.unloading_ratio
            d_history = d_2.flatten()
            d_arr = np.hstack((d_1, d_history))
            d_arr = np.hstack([np.linspace(d_arr[i], d_arr[i + 1], number_of_increments)
                               for i in range(len(d_arr) - 1)])

            return d_arr

        if self.loading_type == "Cyclic" and self.amplitude_type == "Variable_Amplitude" and self.loading_levels == "Flexible_loading_levels(LS3,LS4)":

            d_0 = np.zeros(1)

            d_1 = np.linspace(0, self.max_monotonic *
                              self.S_max_1, self.cycles_n_1 * 2)
            d_1.reshape(-1, 2)[:, 0] = self.max_monotonic * self.S_max_1
            d_1.reshape(-1, 2)[:, 1] = self.max_monotonic * \
                self.S_min_1
            d_history_1 = d_1.flatten()

            d_2 = np.linspace(0, self.max_monotonic *
                              self.S_max_2, self.cycles_n_2 * 2)
            d_2.reshape(-1, 2)[:, 0] = self.max_monotonic * self.S_max_2
            d_2.reshape(-1, 2)[:, 1] = self.max_monotonic *  \
                self.S_min_2
            d_history_2 = d_2.flatten()

            d_3 = np.linspace(0, self.max_monotonic *
                              self.S_max_3, self.cycles_n_3 * 2)
            d_3.reshape(-1, 2)[:, 0] = self.max_monotonic * self.S_max_3
            d_3.reshape(-1, 2)[:, 1] = self.max_monotonic *  \
                self.S_min_3
            d_history_3 = d_3.flatten()

            d_4 = np.linspace(0, self.max_monotonic *
                              self.S_max_4, self.cycles_n_4 * 2)
            d_4.reshape(-1, 2)[:, 0] = self.max_monotonic * self.S_max_4
            d_4.reshape(-1, 2)[:, 1] = self.max_monotonic *  \
                self.S_min_4
            d_history_4 = d_4.flatten()

            d_5 = np.linspace(0, self.max_monotonic *
                              self.S_max_5, self.cycles_n_5 * 2)
            d_5.reshape(-1, 2)[:, 0] = self.max_monotonic * self.S_max_5
            d_5.reshape(-1, 2)[:, 1] = self.max_monotonic *  \
                self.S_min_5
            d_history_5 = d_5.flatten()

            d_arr_1 = np.hstack(
                (d_0, d_history_1, d_history_2, d_history_3, d_history_4, d_history_5))

            d_arr = d_arr_1

            n = self.number_of_repeted_blocks

            if n > 1:
                for i in range(1, n):
                    d_arr = np.hstack(
                        (d_arr, d_history_1, d_history_2, d_history_3, d_history_4, d_history_5))

            d_arr = np.hstack([np.linspace(d_arr[i], d_arr[i + 1], number_of_increments)
                               for i in range(len(d_arr) - 1)])

            self.number_of_cycles = n * \
                (self.cycles_n_1 + self.cycles_n_2 +
                 self.cycles_n_3 + self.cycles_n_4 + self.cycles_n_5)

            return d_arr

        if self.loading_type == "Cyclic" and self.amplitude_type == "Variable_Amplitude" and self.loading_levels == "Increased_blocks(LS2)":

            d_arr_1 = np.zeros(1)

            for i in range(1, self. number_of_levels + 1, 1):

                d_i = np.linspace(0, self.max_monotonic *
                                  (self.s_max_first + (i - 1.0) * (self.s_max_last -
                                                                   self.s_max_first) / (self.number_of_levels - 1.0)),
                                  self.number_of_cycles_each_level * 2)
                d_i.reshape(-1, 2)[:, 0] = self.max_monotonic * (
                    self.s_max_first + (i - 1.0) * (self.s_max_last - self.s_max_first) / (self.number_of_levels - 1.0))
                d_i.reshape(-1, 2)[:, 1] = self.max_monotonic * self.s_min_all

                d_history_i = d_i.flatten()

                d_arr_1 = np.hstack((d_arr_1, d_history_i))

            d_arr = d_arr_1

            n = self.number_of_repeted_blocks
            for i in range(1, n):
                d_arr = np.hstack(
                    (d_arr, d_arr_1[1:]))
            d_arr = np.hstack([np.linspace(d_arr[i], d_arr[i + 1], number_of_increments)
                               for i in range(len(d_arr) - 1)])

            self.number_of_cycles = n * \
                (self.number_of_levels * self.number_of_cycles_each_level)

            return d_arr

    time_func = Property(
        depends_on='loading_type ,amplitude_type, loading_levels ,maximum_loading,\
        S_max_1,  S_max_2, S_max_3, S_max_4,\
        S_max_5, t_max ,d_t, d_array,max_monotonic, number_of_levels,\
        s_max_first ,s_min_all ,s_max_last , number_of_cycles_each_level, number_of_repeted_blocks')

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
        ax.plot(x, self.time_func(x), 'k')
        ax.set_xlabel('time')
        ax.set_ylabel('displacement')
        self.figure.canvas.draw()

    view = View(VGroup(Group(Item('loading_type'),
                             Item('maximum_loading'),
                             Item('number_of_cycles'),
                             Item('amplitude_type'),
                             Item('loading_levels'),
                             Item('unloading_ratio'), show_border=True, label='Cyclic load inputs'),

                       Group(Item('d_t'),
                             Item('t_max'),
                             show_border=True,
                             label='Solver Settings'),

                       Group(Item('max_monotonic'),
                             HGroup(Item('number_of_levels'),
                                    Item('number_of_cycles_each_level')),
                             HGroup(Item('s_max_first'),
                                    Item('s_max_last')),
                             Item('s_min_all'), show_border=True, label='Increased Loading (LS2)'),

                       Group(HGroup(Item('cycles_n_1'),
                                    Item('S_max_1')),
                             Item('S_min_1'),
                             HGroup(
                           Item('cycles_n_2'),
                           Item('S_max_2')), Item('S_min_2'),
        HGroup(Item('cycles_n_3'), Item('S_max_3')),
        Item('S_min_3'),
        HGroup(Item('cycles_n_4'), Item('S_max_4')),
        Item('S_min_4'),
        HGroup(Item('cycles_n_5'), Item('S_max_5')),
        Item('S_min_5'),
        Item('number_of_repeted_blocks'),
        show_border=True, label='flexible Loading Blocks (LS3,LS4)')),



        Group(Item('update', label='Plot Loading scenario')),
        Item('figure', editor=MPLFigureEditor(),
             dock='horizontal', show_label=False),
        Item('time', label='t/T_max'))


class AllicheConcreteFatigueModel(BMCSTreeNode):

    node_name = Str('Concrete_Fatigue_Alliche_model')

    tree_node_list = List([])

    def _tree_node_list_default(self):
        return [self.material, self.loading_scenario]

    material = Instance(Material)

    def _material_default(self):
        return Material()

    loading_scenario = Instance(LoadingScenario)

    def _loading_scenario_default(self):
        return LoadingScenario()

    mats_eval = Instance(MATSEvalAllicheFatigue)

    figure = Instance(Figure)

    def _figure_default(self):
        figure = Figure(facecolor='white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    def plot(self, figure):
        self.mats_eval.lamda = self.material.lamda
        self.mats_eval.mu = self.material.mu
        self.mats_eval.alpha = self.material.alpha
        self.mats_eval.beta = self.material.beta
        self.mats_eval.g = self.material.g
        self.mats_eval.C0 = self.material.C0
        self.mats_eval.C1 = self.material.C1
        self.mats_eval.K = self.material.K
        self.mats_eval.n = self.material.n

        s_arr = self.loading_scenario._get_d_array()
        t_arr = np.linspace(0, 1, len(s_arr))
        m = self.loading_scenario.t_max / self.loading_scenario.d_t
        n = len(s_arr)

        sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, inc = self.mats_eval.get_stress_strain(
            s_arr * -1.0)

        # Fig1 - Loading history
        ax1 = figure.add_subplot(221)
        ax1.cla()
        ax1.plot(
            t_arr[0:inc], abs(sigma_arr[0:inc]), 'k', linewidth=1, alpha=1.0)
        ax1.set_title('loading history')
        ax1.set_xlabel('Time')
        ax1.set_ylabel('$\sigma_{1}$')

        # Fig2 - stress-strain curve
        ax2 = figure.add_subplot(222)
        ax2.plot(abs(eps_1_arr[0:inc]), abs(
            sigma_arr[0:inc]), 'k', linewidth=1, alpha=1.0)
        ax2.set_title('$ \epsilon_{11} - \sigma_{11}$')
        ax2.set_xlabel('$\epsilon_{11}$')
        ax2.set_ylabel('$\sigma_{11}$[MPa]')

        # Fig3 - damage evolution
        ax3 = figure.add_subplot(223)

        w = np.zeros(n)
        cycle = np.zeros(n)
        for i in range(0, n, 1):
            idx = m + 2 * i * m - 1
            if idx <= len(w_arr[0:inc]):
                idx = idx
            else:
                idx = m + 2 * (i - 1.0) * m - 1
                break

            w[i] = w_arr[int(idx)]
            cycle[i] = i + 1

        ax3.plot(cycle[0:i], w[0:i], 'k', linewidth=1, alpha=1)
        ax3.set_xlabel('number of cycles')
        ax3.set_ylabel('Damage')

        # Fig3 - fatigue creep curve
        ax4 = figure.add_subplot(224)
        eps_1_max = np.zeros(n)
        eps_1_min = np.zeros(n)
        cycle = np.zeros(n)
        for i in range(0, n, 1):
            idx_1 = m + 2 * i * m - 1
            idx_2 = 2 * i * m
            if idx_1 <= len(eps_1_arr[0:inc]):
                idx_1 = idx_1
            else:
                idx_1 = m + 2 * (i - 1.0) * m - 1
                break

            if idx_2 <= len(eps_1_arr[0:inc]):
                idx_2 = idx_2
            else:
                idx_2 = 1 * (i - 1.0) * m
                break

            eps_1_max[i] = eps_1_arr[int(idx_1)]
            eps_1_min[i] = eps_1_arr[int(idx_2)]
            cycle[i] = i + 1

        ax4.plot(cycle[0:i], abs(eps_1_max[0:i]), 'k', linewidth=1, alpha=1)
        ax4.plot(cycle[1:i], abs(eps_1_min[1:i]), 'k', linewidth=1, alpha=1)

        # experimental results IMB
        n_max_1 = np.loadtxt(
            r'E:\Publishing\FIB_2018_(1)\results\Alliche\calibration\exp\LS2_Nmax_1.txt')
        n_min_1 = np.loadtxt(
            r'E:\Publishing\FIB_2018_(1)\results\Alliche\calibration\exp\LS2_Nmin_1.txt')
        s_max_1 = np.loadtxt(
            r'E:\Publishing\FIB_2018_(1)\results\Alliche\calibration\exp\LS2_Smax_1.txt')
        s_min_1 = np.loadtxt(
            r'E:\Publishing\FIB_2018_(1)\results\Alliche\calibration\exp\LS2_Smin_1.txt')

        n_max_2 = np.loadtxt(
            r'E:\Publishing\FIB_2018_(1)\results\Alliche\calibration\exp\LS2_Nmax_2.txt')
        n_min_2 = np.loadtxt(
            r'E:\Publishing\FIB_2018_(1)\results\Alliche\calibration\exp\LS2_Nmin_2.txt')
        s_max_2 = np.loadtxt(
            r'E:\Publishing\FIB_2018_(1)\results\Alliche\calibration\exp\LS2_Smax_2.txt')
        s_min_2 = np.loadtxt(
            r'E:\Publishing\FIB_2018_(1)\results\Alliche\calibration\exp\LS2_Smin_2.txt')

        ax4.plot(n_max_1[1:], s_max_1[1:] / 300, ":k", label='Smax')
        ax4.plot(n_min_1[1:], s_min_1[1:] / 300, ":k", label='Smin')
        ax4.plot(n_max_2[1:], s_max_2[1:] / 300, ":k", label='Smax')
        ax4.plot(n_min_2[1:], s_min_2[1:] / 300, ":k", label='Smin')

        ax4.set_xlabel('number of cycles')
        ax4.set_ylabel('max $\epsilon_{11}$')

        #======================================================================
        # saving results
        #======================================================================
        #eps_record_1 = np.zeros(1)
        #sigma_record = np.zeros(1)
        #eps_record_2 = np.zeros(1)
        N_record = np.zeros(1)
        w_record = np.zeros(1)
        eps_max_record = np.zeros(1)
        eps_min_record = np.zeros(1)
        #stiffness_record = np.zeros(1)

        for i in range(0, i, 1):
            #eps_record_1 = np.vstack((eps_record_1, eps_1_arr[i]))
            #sigma_record = np.vstack((sigma_record, sigma_arr[i]))
            #eps_record_2 = np.vstack((eps_record_2, eps_2[i]))
            N_record = np.vstack((N_record, cycle[i]))
            eps_max_record = np.vstack((eps_max_record, eps_1_max[i]))
            eps_min_record = np.vstack((eps_min_record, eps_1_min[i]))
            w_record = np.vstack((w_record, w[i]))
    #         stiffness_record = np.vstack(
    #             (stiffness_record, maximum_stress / eps_1[i]))

        np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\eps_1.txt',
                   np.transpose(eps_1_arr), delimiter=" ", fmt="%s")
        np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\sigma_1.txt',
                   np.transpose(sigma_arr), delimiter=" ", fmt="%s")
        np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\eps_max.txt',
                   np.transpose(eps_max_record), delimiter=" ", fmt="%s")
        np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\eps_min.txt',
                   np.transpose(eps_min_record), delimiter=" ", fmt="%s")
        np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\N.txt',
                   np.transpose(N_record), delimiter=" ", fmt="%s")
        np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\w.txt',
                   np.transpose(w_record), delimiter=" ", fmt="%s")
    #     np.savetxt(r'E:\Models_Implementation\Concrete Fatigue models\Alliche_2004\Results\loading_sequence\creep_fatigue\stiffness.txt',
    #                np.transpose(stiffness_record), delimiter=" ", fmt="%s")


if __name__ == '__main__':

    material_node = TreeNode(node_for=[Material],
                             auto_open=False,
                             children='tree_node_list',
                             label='node_name',
                             view='tree_view',
                             )

    loading_scenario_node = TreeNode(node_for=[LoadingScenario],
                                     auto_open=True,
                                     children='tree_node_list',
                                     label='node_name',
                                     view='tree_view',
                                     )

    bond_slip_model_node = TreeNode(node_for=[AllicheConcreteFatigueModel],
                                    auto_open=True,
                                    children='tree_node_list',
                                    label='node_name',
                                    view='tree_view',
                                    #menu=Menu(plot_self, NewAction),
                                    )

# =========================================================================
# List of all custom nodes
# =========================================================================

    custom_node_list = [material_node, loading_scenario_node,
                        bond_slip_model_node]

    model = AllicheConcreteFatigueModel(mats_eval=MATSEvalAllicheFatigue())
    w = BMCSWindow(root=model)
    w.configure_traits()
