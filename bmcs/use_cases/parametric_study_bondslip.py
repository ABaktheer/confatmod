'''
Created on 02.01.2017

@author: abaktheer
'''
from traits.api import \
    Instance, Property, \
    List, Str, Trait, Button
from traitsui.api import \
    View, Item, UItem, VGroup, HGroup, spring
from view.ui.bmcs_tree_node import \
    BMCSTreeNode
from view.window.bmcs_window import \
    BMCSWindow

from matmod.bond_slip_model_fatigue import \
    BondSlipModel, Material, LoadingScenario
from matmod.mats_bondslip_fatigue import \
    MATSEvalFatigue
import matplotlib.gridspec as gridspec
from utils.keyref import \
    KeyRef


class UCPStudyElement(BMCSTreeNode):
    '''Class controlling plotting options
    for an instance
    '''
    node_name = Str('<unnamed>')

    color = Trait('black', dict(black='k',
                                cyan='c',
                                green='g',
                                blue='b',
                                yellow='y',
                                magneta='m',
                                red='r',
                                gray_1='0.35',
                                gray_2='0.5',
                                gray_3='0.65')
                  )

    linestyle = Trait('solid', dict(solid='-',
                                    dashed='--',
                                    dash_dot='-.',
                                    dotted=':')
                      )

    linewidth = Trait('medium', dict(thin=1.0,
                                     medium=1.5,
                                     thik=2.0,
                                     )
                      )

    alpha = Trait('full', dict(full=1.0,
                               half=0.5
                               )
                  )

    tree_view = View(VGroup(Item('node_name', label='label'),
                            Item('linestyle'),
                            Item('linewidth'),
                            Item('alpha'),
                            Item('color'),
                            label='Plotting options'))

    def plot(self, fig):
        self.content.plot(fig, color=self.color_, linestyle=self.linestyle_, linewidth=self.linewidth_, alpha=self.alpha_,
                          label=self.node_name)

    def plot_ax(self, ax1, ax2, ax3, ax4, ax5, ax6):
        self.content.plot_custom(ax1=ax1, ax2=ax2, ax3=ax3, ax4=ax4, ax5=ax5, ax6=ax6, color=self.color_, linestyle=self.linestyle_, linewidth=self.linewidth_, alpha=self.alpha_,
                                 label=self.node_name)


class UCPStudyElementBMCS(UCPStudyElement):
    node_name = '<unnamed bond_slip>'

    tree_node_list = List(Instance(BMCSTreeNode))

    def _tree_node_list_default(self):
        return [BondSlipModel(mats_eval=MATSEvalFatigue())]

    content = Property(depends_on='tree_node_list')

    def _get_content(self):
        return self.tree_node_list[0]

    def _set_content(self, val):
        self.tree_node_list = [val]


class UCParametricStudy(BMCSTreeNode):
    node_name = Str('Parametric study')

    element_to_add = Trait(
        'BondSlipModel', {'BondSlipModel':   UCPStudyElementBMCS})

    add_element = Button('Add')

    def _add_element_fired(self):
        self.append_node(self.element_to_add_())

    tree_view = View(HGroup(UItem('element_to_add', springy=True),
                            UItem('add_element')),
                     spring
                     )

    tree_node_list = List(Instance(BMCSTreeNode))

    def _tree_node_list_default(self):
        return []

    def plot(self, fig):
        ax1 = fig.add_subplot(231)
        ax2 = fig.add_subplot(232)
        ax3 = fig.add_subplot(233)
        ax4 = fig.add_subplot(234)
        ax5 = fig.add_subplot(235)
        ax6 = fig.add_subplot(236)

        for node in self.tree_node_list:
            node.plot_ax(ax1, ax2, ax3, ax4, ax5, ax6)


bond_slip_ps = UCParametricStudy()
bond_slip_ps.element_to_add = 'BondSlipModel'
bond_slip_ps.add_element = True
bond_slip_ps.add_element = True

ucc = BMCSTreeNode()
ucc.tree_node_list.append(bond_slip_ps)

mxn_ps_view = BMCSWindow(root=ucc)
mxn_ps_view.configure_traits()
