from neuron import h, gui
from neuron.units import ms, mV
import matplotlib.pyplot as plt
h.load_file('stdrun.hoc')

class Cell:
    def __init__(self, gid, x, y, z, theta):
        self._gid = gid
        self._setup_morphology()
        self.all = self.soma.wholetree()
        self._setup_biophysics()
        self.x = self.y = self.z = 0
        h.define_shape()
        self._rotate_z(theta)
        self._set_position(x, y, z)

        self._spike_detector = h.NetCon(self.soma(0.5)._ref_v, None, sec = self.soma)
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

class Sensory(Cell):
    name = "Sensory"
    def _setup_morphology(self):
        self.soma = h.Section(name = 'soma', cell = self)
        self.peripheral_axon = h.Section(name = 'peripheral_axon', cell = self)
        self.central_axon = h.Section(name = 'central_axon', cell = self)
        self.peripheral_axon.connect(self.soma(0))
        self.central_axon.connect(self.soma(0))
        self.soma.L = self.soma.diam = 25

        self.peripheral_axon.L = 5000
        self.peripheral_axon.diam = 1.0
        self.peripheral_axon.nseg = 101

        self.central_axon.L = 500
        self.central_axon.diam = 1.0
        self.central_axon.nseg = 21
    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100
            sec.cm = 1
        self.soma.insert('pas')
        self.central_axon.insert('hh')
        self.peripheral_axon.insert('hh')
        self.stim = h.IClamp(self.peripheral_axon(1))
        self.stim.delay = 5     # ms
        self.stim.dur = 1       # ms
        self.stim.amp = 0.2     # nA

        self.t = h.Vector().record(h._ref_t)
        self.v_soma = h.Vector().record(self.soma(0.5)._ref_v)
        self.v_peripheral = h.Vector().record(self.peripheral_axon(0.5)._ref_v)
        self.v_central = h.Vector().record(self.central_axon(0.5)._ref_v)
    def set_stim(self, delay=5, dur=1, amp=0.3):
        self.stim.delay = delay
        self.stim.dur = dur
        self.stim.amp = amp


sensory = Sensory(1,0,0,0,0)
sensory.set_stim(delay=2, dur=1, amp=0.3)

h.finitialize(-65)
h.continuerun(40)

import matplotlib.pyplot as plt
plt.plot(sensory.t, sensory.v_peripheral, label='Peripheral Axon')
plt.plot(sensory.t, sensory.v_soma, label='Soma')
plt.plot(sensory.t, sensory.v_central, label='Central Axon')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
plt.legend()
plt.title("Sensory Neuron Firing")
plt.grid()
plt.show()
