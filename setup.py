from neuron import h, gui
from neuron.units import ms, mV
import matplotlib.pyplot as plt
h.load_file('stdrun.hoc')

# --------------------------------- change variables here --------------------------------------------------------------

V_REST = -65  

spike_detector_loc = 0.5


# Node parameters
node_cm = 1.0
node_Ra = 100
node_L = 1
node_nseg = 1


# Internode parameters
#internode_cm = 0.04
internode_cm = 0.01
internode_g_pas = 1e-5
internode_e_pas = V_REST
internode_Ra = 100
internode_nseg = 31

# Sensory neuron geometry
p_L = 800              # Peripheral axon total length
p_int_L = 500           # Peripheral internode length
p_diam = 14            # Peripheral axon diameter
c_L = 500               # Central axon total length
c_int_L = 90            # Central internode length
c_diam = 1.0            # Central axon diameter
sensory_soma_L_diam = 25

# Channel properties (HH mechanism)
gnabar = 0.3          # Sodium conductance
gkbar = 0.1           # Potassium conductance
gl = 0.0003            # Leak conductance
el = V_REST              # Leak reversal potential

# Passive dendrites
g_pas_dend = 1e-4
e_pas_dend = V_REST

# Soma/Dend parameters (motor neuron)
motor_soma_L_diam = 30
motor_dend_L = 200
motor_dend_diam = 2

# Interneuron parameters
interneuron_soma_L_diam = 20
interneuron_dend_L = 100
interneuron_dend_diam = 2
interneuron_axon_L = 100
interneuron_axon_int_L = 15
interneuron_axon_diam = 4



# ------------------------------------------- Class definitions --------------------------------------------------------
class Cell:
    
    def __init__(self, gid, x, y, z, theta):
        self.active_sections = []
        self._gid = gid
        self._setup_morphology()
        self.all = self.soma.wholetree()
        self._setup_biophysics()
        self.x = self.y = self.z = 0
        h.define_shape()
        self._rotate_z(theta)
        self._set_position(x, y, z)
        

        self._spike_detector = h.NetCon(self.soma(spike_detector_loc)._ref_v, None, sec = self.soma)
        self.spike_times = h.Vector()
        self._spike_detector.record(self.spike_times)
        self._ncs = []

    def __repr__(self):
        return '{}[{}]'.format(self.name, self._gid)
    
    def _set_position(self, x, y, z):
        for sec in self.all:
            for i in range(sec.n3d()):
                sec.pt3dchange(i, x-self.x+sec.x3d(i), y-self.y+sec.y3d(i), z-self.z+sec.z3d(i), sec.diam3d(i))
        self.x, self.y, self.z = x, y, z
    def _rotate_z(self, theta):
        for sec in self.all:
            for i in range(sec.n3d()):
                x = sec.x3d(i)
                y = sec.y3d(i)
                c = h.cos(theta)
                s = h.sin(theta)
                xprime = x*c - y*s
                yprime = x*s + y*c
                sec.pt3dchange(i, xprime, yprime, sec.z3d(i), sec.diam3d(i))
    def myelinated(self, L_axon, internode_L, diam, nnodes=-1):
        axon_sections = []
        total_length = 0
        n_nodes = 0
        i = 0
        while total_length < L_axon and (n_nodes < nnodes or nnodes == -1):
            node = h.Section(name=f'node_{i}')
            self.active_sections.append(node)
            node.insert('hh')
            node.cm = node_cm
            node.Ra = node_Ra
            node.L = node_L
            total_length += node.L
            node.diam = diam
            node.nseg = node_nseg
            if (nnodes != -1):
                n_nodes += 1
            
            internode = h.Section(name=f'internode_{i}')
            internode.insert('pas')
            internode.cm = internode_cm
            internode.g_pas = internode_g_pas
            internode.e_pas = internode_e_pas
            internode.Ra = internode_Ra
            internode.L = internode_L
            internode.diam = diam
            internode.nseg = internode_nseg
            total_length += internode.L

            if axon_sections:
                internode.connect(axon_sections[-1](1))
            node.connect(internode(1))
            axon_sections.append(internode)
            axon_sections.append(node)
            i+= 1
        return axon_sections

class Sensory(Cell):
    name = "Sensory"
    
    def _setup_morphology(self):
        self.soma = h.Section(name = 'soma', cell = self)
        
       #self.peripheral_axon = self.myelinated(p_L, p_int_L, p_diam)
        self.central_axon = self.myelinated(c_L, c_int_L, c_diam)
        #self.peripheral_axon[0].connect(self.soma(0))
        self.central_axon[0].connect(self.soma(1))
        self.soma.L = self.soma.diam = sensory_soma_L_diam

    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100
            sec.cm = 1
        self.soma.insert('hh')
        """for sec in self.active_sections + [self.soma]:
            for seg in sec:
                seg.hh.gnabar = gnabar
                seg.hh.gkbar = gkbar
                seg.hh.gl = gl
                seg.hh.el = el"""
        self.stimI = h.IClamp(self.soma(0.5))
        #for sec in reversed(self.peripheral_axon):
        #    if sec in self.active_sections:
        #        self.stimI = h.IClamp(sec(0.5))
        #        break

        #self.stim = h.NetStim()
        #self.stim.delay = 5     # ms
        #self.stim.dur = 5       # ms
        #self.stim.amp = 2     # nA

        self.t = h.Vector().record(h._ref_t)
        self.v_soma = h.Vector().record(self.soma(0.5)._ref_v)
        #self.v_peripheral = h.Vector().record(self.peripheral_axon[-1](0.5)._ref_v)

        self.v_central = h.Vector().record(self.central_axon[-1](0.5)._ref_v)
    def set_stim(self, delay=5, dur=1, amp=0.3):
        self.stimI.delay = delay
        self.stimI.dur = dur
        self.stimI.amp = amp

class Motor(Cell):
    name = "Motor"
    def _setup_morphology(self):
        self.soma = h.Section(name = 'soma', cell = self)
        self.soma.L = self.soma.diam = 30
        self.dend = h.Section(name="dend", cell=self)
        self.dend.L = motor_dend_L
        self.dend.diam = motor_dend_diam
        self.dend.connect(self.soma(0))
        
        #self.axon = self.myelinated(L_axon=5100, internode_L = 1000, diam=10)
        #self.axon[0].connect(self.soma(1))
        self.axon = self.myelinated(L_axon=5000, internode_L=1000, diam=1.0)
        self.axon[0].connect(self.soma(1))


    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100
            sec.cm = 1
        self.soma.insert("hh")
        self.dend.insert('hh')
        #self.soma.insert('kad')
        #self.soma.gkbar_kad = 0.005
        for seg in self.dend:
            seg.hh.gnabar = 0.04
            seg.hh.gkbar  = 0.01
            seg.hh.gl     = gl
            seg.hh.el     = el

        for sec in self.active_sections + [self.soma]:
            for seg in sec:
                seg.hh.gnabar = gnabar
                seg.hh.gkbar = gkbar
                seg.hh.gl = gl
                seg.hh.el = el

        #self.stim = h.IClamp(self.dend(0.5))
        #self.stim.delay = 5     # ms
        #self.stim.dur = 1       # ms
        #self.stim.amp = 0.05     # nA

        self.t = h.Vector().record(h._ref_t)
        self.v_soma = h.Vector().record(self.soma(0.5)._ref_v)
        self.v_dend = h.Vector().record(self.dend(0.5)._ref_v)
        ind = -1
        #self.v_axon = h.Vector().record(self.axon[ind](0.5)._ref_v)

    def set_stim(self, delay=5, dur=1, amp=0.3):
        self.stim.delay = delay
        self.stim.dur = dur
        self.stim.amp = amp


class MyelinatedInterneuron(Cell):
    name = "MyelinatedInterneuron"

    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.soma.L = self.soma.diam = 20

        self.dendrite = h.Section(name='dendrite', cell=self)
        self.dendrite.L = 100
        self.dendrite.diam = 2
        self.dendrite.connect(self.soma(0))

        self.axon = self.myelinated(L_axon=interneuron_axon_L, internode_L=interneuron_axon_int_L, diam=interneuron_axon_diam)
        self.axon[0].connect(self.soma(1))  # connect first internode to soma

        self.active_sections.extend([self.dendrite])  # axon parts already in active_sections

    def _setup_biophysics(self):
        for sec in [self.soma, self.dendrite]:
            sec.Ra = 100
            sec.cm = 1
            sec.insert('hh')
            
            for seg in sec:
                seg.hh.gnabar = gnabar
                seg.hh.gkbar = gkbar
                seg.hh.gl = gl
                seg.hh.el = el

        #self.soma.insert('kad')
        #self.soma.gkbar_kad = 0.01
        self.t = h.Vector().record(h._ref_t)
        self.v_soma = h.Vector().record(self.soma(0.5)._ref_v)
        self.v_dend = h.Vector().record(self.dendrite(0.5)._ref_v)
        self.v_axon = h.Vector().record(self.axon[3](0.5)._ref_v)  # record mid-axon

        # Stimulus at dendrite
        #self.stim = h.IClamp(self.dendrite(0.5))
        #self.stim.delay = 5
        #self.stim.dur = 5
        #self.stim.amp = 0.5

    def set_stim(self, delay=5, dur=1, amp=0.5):
        self.stim.delay = delay
        self.stim.dur = dur
        self.stim.amp = amp


