from neuron import h, gui
from neuron.units import ms, mV
import matplotlib.pyplot as plt
h.load_file('stdrun.hoc')
from setup import Cell, Sensory, Motor, MyelinatedInterneuron