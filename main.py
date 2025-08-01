
#TO DO
#make sure inhibition/excitation works as necessary
#tweak sensory input so that it is more realistic
#slow potassium channels on modeldb

from neuron import h, gui
from neuron.units import ms, mV
import matplotlib.pyplot as plt
h.load_file('stdrun.hoc')
#h.nrn_load_dll("x86_64/.libs/libnrnmech.so")
from setup import Cell, Sensory, Motor, MyelinatedInterneuron

V_REST = -65  

spike_detector_loc = 0.5


# Synaptic parameters
syn_si_tau1 = 0.3 #
syn_si_tau2 = 2.0 #
syn_si_e = 0
syn_si_threshold = -20
syn_si_delay = 1 #
syn_si_weight = 0.03 #

syn_im_tau1 = 0.2
syn_im_tau2 = 5.0
syn_im_e = -80
syn_im_threshold = -20
syn_im_delay = 1.0
syn_im_weight = 0.015

syn_sm_tau1 = 1.5
syn_sm_tau2 = 2.0
syn_sm_e = 0
syn_sm_threshold = -20
syn_sm_delay = 1.18
syn_sm_weight = 0.02

# ====================================== Simulation ------------------------------------------------------------------------

sensory = Sensory(0,0,0,0,0)
sensory.set_stim(delay=2, dur=100, amp=1)

"""
sensory.stim.start = 2    # ms
sensory.stim.number = 10
sensory.stim.interval = 1  # not used if number = 1
sensory.stim.noise = 0.2

syn = h.ExpSyn(sensory.soma(0.5))
syn.tau = 0.5
syn.e = 0  # Excitatory

nc = h.NetCon(sensory.stim, syn)
nc.delay = 0
nc.weight[0] = 1  # Adjust strength
"""

        
interneuron = MyelinatedInterneuron(3, 50, 0, 0, 0)
#interneuron.set_stim(delay=2, dur=5, amp=1)

motor = Motor(0, 0, 0, 0, 0)
#motor.stim.amp = 3  # Stronger if no spike

syn_si = h.Exp2Syn(interneuron.dendrite(0.5))
syn_si.tau1 = syn_si_tau1
syn_si.tau2 = syn_si_tau2
syn_si.e = syn_si_e

nc_si = h.NetCon(sensory.soma(spike_detector_loc)._ref_v, syn_si, sec=sensory.soma)
nc_si.threshold = syn_si_threshold
nc_si.delay = syn_si_delay
nc_si.weight[0] = syn_si_weight


#right now modeling inhibition
syn_im = h.Exp2Syn(motor.dend(0.75))
syn_im.tau1 = syn_im_tau1
syn_im.tau2 = syn_im_tau2
syn_im.e = syn_im_e
nc_im = h.NetCon(interneuron.soma(spike_detector_loc)._ref_v, syn_im, sec = interneuron.soma)
nc_im.threshold = syn_im_threshold
nc_im.delay = syn_im_delay
nc_im.weight[0] = syn_im_weight


syn_sm = h.Exp2Syn(motor.dend(0.3))  # target more proximal part of dendrite
syn_sm.tau1 = syn_sm_tau1
syn_sm.tau2 = syn_sm_tau2
syn_sm.e = syn_sm_e

nc_sm = h.NetCon(sensory.soma(spike_detector_loc)._ref_v, syn_sm, sec=sensory.soma)
nc_sm.threshold = syn_sm_threshold
nc_sm.delay = syn_sm_delay
nc_sm.weight[0] = syn_sm_weight

stim_amp = h.Vector().record(sensory.stimI._ref_i)

h.finitialize(V_REST)
h.continuerun(100)

# ----------------------------------------------- Plots ------------------------------------------------------------------
plt.subplot(4,1,1)
#plt.plot(sensory.t, sensory.v_peripheral, label='Peripheral Axon')
plt.plot(sensory.t, sensory.v_soma, label='Soma')
plt.plot(sensory.t, sensory.v_central, label='Central Axon')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
plt.legend()
plt.title("Sensory Neuron Firing")
plt.grid()

plt.subplot(4,1,2)
plt.plot(interneuron.t, interneuron.v_soma, label='Soma')
plt.plot(interneuron.t, interneuron.v_dend, label='Dendrite') 
plt.plot(interneuron.t, interneuron.v_axon, label='Myelinated Axon')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
plt.legend()
plt.title("Myelinated Interneuron Firing")
plt.grid()


plt.subplot(4,1,3)
#plt.plot(motor.t, motor.v_dend, label="Dendrite")
plt.plot(motor.t, motor.v_soma, label="Soma")
#plt.plot(motor.t, motor.v_axon, label="Axon terminal")
plt.legend()
plt.xlabel("Time (ms)")
plt.ylabel("Membrane Potential (mV)")
plt.title("Motor Neuron Response")
plt.grid()


plt.subplot(4,1,4)
plt.plot(sensory.t, stim_amp, label="Injected Current")
plt.legend()
plt.xlabel("Time (ms)")
plt.ylabel("Injected Current (nA)")
plt.title("Injected Stimulus")
plt.grid()
plt.show()



print(list(sensory.spike_times))
