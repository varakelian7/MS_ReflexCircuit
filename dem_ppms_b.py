from neuron import h, gui
from neuron.units import ms, mV
import matplotlib.pyplot as plt
import numpy as np
h.load_file('stdrun.hoc')
from setup import Sensory, Motor, MyelinatedInterneuron
h.nrn_load_dll('./arm64/.libs/libnrnmech.dylib')
import random as rd


V_REST = -65  

DEMYE_PCT = 21

SIM_DUR = 2000

spike_detector_loc = 0.5

# Synaptic parameters
syn_si_tau1 = 0.3 #
syn_si_tau2 = 2.0 #
syn_si_e = 0
syn_si_threshold = -20
syn_si_delay = 0 #
syn_si_weight = 0.03 #0.03

syn_im_tau1 = 0.2
syn_im_tau2 = 5.0
syn_im_e = -80
syn_im_threshold = -20
syn_im_delay = 0 #1
syn_im_weight = 0.02 #0.015

syn_sm_tau1 = 1.5
syn_sm_tau2 = 2.0
syn_sm_e = 0
syn_sm_threshold = -20
syn_sm_delay = 0
syn_sm_weight = 0.02 #0.02

# ====================================== Simulation ------------------------------------------------------------------------

sensory = Sensory(0,0,0,0,0)


sensory.set_stim(delay=2, dur=SIM_DUR, amp=.7)

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

#Demyelination changes:
s_internodes = []
for internode in sensory.central_axon:
    if internode not in sensory.active_sections:
        s_internodes.append(internode)
m_internodes = []
for internode in motor.axon:
    if internode not in motor.active_sections:
        m_internodes.append(internode)
i_internodes = []
for internode in interneuron.axon:
    if internode not in interneuron.active_sections:
        i_internodes.append(internode)

for section in s_internodes:
    section.insert('hh')  # Add Hodgkin-Huxley channels
    section.gnabar_hh = 0.012
    section.gkbar_hh = 0.0036
    section.cm = 0.3
    section.Ra = 50*1.25
    section.g_pas = 1e-5
    i = 0
    for seg in section:
        if i % round(100/DEMYE_PCT) == 0:  # 1 in every 100/DEMYE_PCT should get it
            seg.gnabar_hh = 0.12
            seg.gkbar_hh = 0.036
            seg.cm = 0.3*1.25
            seg.g_pas = 0.001
        i+=1


for section in i_internodes:
    section.insert('hh')  # Add Hodgkin-Huxley channels
    section.gnabar_hh = 0.012
    section.gkbar_hh = 0.0036
    section.cm = 0.3
    section.Ra = 50*1.25
    section.g_pas = 1e-5
    i = 0
    for seg in section:
        if i % round(100/DEMYE_PCT) == 0:  # 1 in every 100/DEMYE_PCT should get it
            seg.gnabar_hh = 0.12
            seg.gkbar_hh = 0.036
            seg.cm = 0.3*1.25
            seg.g_pas = 0.001
        i+=1

for section in m_internodes:
    section.insert('hh')  # Add Hodgkin-Huxley channels
    section.gnabar_hh = 0.012
    section.gkbar_hh = 0.0036
    section.cm = 0.3
    section.Ra = 50*1.25
    section.g_pas = 1e-5
    i = 0
    for seg in section:
        if i % round(100/DEMYE_PCT) == 0:  # 1 in every 100/DEMYE_PCT should get it
            seg.gnabar_hh = 0.12
            seg.gkbar_hh = 0.036
            seg.cm = 0.3*1.25
            seg.g_pas = 0.001
        i+=1


syn_si = h.Exp2Syn(interneuron.dendrite(0.5))
syn_si.tau1 = syn_si_tau1
syn_si.tau2 = syn_si_tau2
syn_si.e = syn_si_e

syn_nmda_si = h.DetAMPANMDA(interneuron.dendrite(0.5))
nc_nmda_si = h.NetCon(sensory.central_axon[-1](0.5)._ref_v, syn_nmda_si, sec=sensory.central_axon[-1])
nc_nmda_si.threshold = syn_si_threshold  
nc_nmda_si.weight[0] = 0.01*1.1



nc_si = h.NetCon(sensory.central_axon[-1](0.5)._ref_v, syn_si, sec=sensory.central_axon[-1])
nc_si.threshold = syn_si_threshold
nc_si.delay = syn_si_delay
nc_si.weight[0] = syn_si_weight


#right now modeling inhibition

syn_im = h.Exp2Syn(motor.dend(0.5))
syn_im.tau1 = syn_im_tau1
syn_im.tau2 = syn_im_tau2
syn_im.e = syn_im_e
nc_im = h.NetCon(interneuron.axon[-1](0.5)._ref_v, syn_im, sec = interneuron.axon[-1])
nc_im.threshold = syn_im_threshold
nc_im.delay = syn_im_delay
nc_im.weight[0] = syn_im_weight

gabab_syn = h.GABAb_S(motor.dend(0.5))  # adjust location as needed
gabab_syn.gmax = 0.01225*1.8  # adjust conductance

# Set up NetCon
gabab_nc = h.NetCon(interneuron.axon[-1](0.5)._ref_v, None, sec=interneuron.axon[-1])
gabab_nc.threshold = syn_im_threshold
gabab_nc.delay = 1
gabab_nc.weight[0] = 0.1*0.95  # this weight will toggle between 0/1 as event flag

# Link the pointer to the weight[0]
h.setpointer(gabab_nc._ref_weight[0], 'pre', gabab_syn)


syn_sm = h.Exp2Syn(motor.dend(0.5))
syn_sm.tau1 = syn_sm_tau1
syn_sm.tau2 = syn_sm_tau2
syn_sm.e = syn_sm_e

syn_nmda_sm = h.DetAMPANMDA(motor.dend(0.5))
nc_nmda_sm = h.NetCon(sensory.central_axon[-1](0.5)._ref_v, syn_nmda_sm, sec=sensory.central_axon[-1])
nc_nmda_sm.threshold = syn_sm_threshold
nc_nmda_sm.weight[0] = 0.01*1.1



nc_sm = h.NetCon(sensory.central_axon[-1](0.5)._ref_v, syn_sm, sec=sensory.central_axon[-1])
nc_sm.threshold = syn_sm_threshold
nc_sm.delay = syn_sm_delay
nc_sm.weight[0] = syn_sm_weight

stim_amp = h.Vector().record(sensory.stimI._ref_i)

#random noise:

s_noise = h.IClamp(sensory.soma(0.5))
s_noise.delay = 0
s_noise.dur = SIM_DUR

i_noise = h.IClamp(interneuron.soma(0.5))
i_noise.delay = 0
i_noise.dur = SIM_DUR

m_noise = h.IClamp(motor.soma(0.5))
m_noise.delay = 0
m_noise.dur = SIM_DUR

dt = 0.1  # ms, same as h.dt
tstop = SIM_DUR  # ms
npts = int(tstop / dt)

# Generate Gaussian noise current: mean = 0, std = 0.05 nA
noise_current = np.random.normal(loc=0.0, scale=0.05, size=npts)
time_vector = np.arange(0, tstop, dt)

vec_i = h.Vector(noise_current)
vec_t = h.Vector(time_vector)

vec_i.play(s_noise._ref_amp, vec_t, 1)  # 1 = continuous interpolation
vec_i.play(i_noise._ref_amp, vec_t, 1)  # 1 = continuous interpolation
vec_i.play(m_noise._ref_amp, vec_t, 1)  # 1 = continuous interpolation


h.finitialize(V_REST)
h.continuerun(SIM_DUR)

# ----------------------------------------------- Plots ------------------------------------------------------------------
plt.subplot(4,1,1)
#plt.plot(sensory.t, sensory.v_peripheral, label='Peripheral Axon')
plt.plot(sensory.t, sensory.v_soma, label='Soma')
#plt.plot(sensory.t, sensory.v_central, label='Central Axon')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
#plt.legend()
plt.title("Sensory Neuron (PPMS)")
plt.grid()

plt.subplot(4,1,2)
plt.plot(interneuron.t, interneuron.v_soma, label='Soma')
#plt.plot(interneuron.t, interneuron.v_dend, label='Dendrite') 
#plt.plot(interneuron.t, interneuron.v_axon, label='Myelinated Axon')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
#plt.legend()
plt.title("Interneuron")
plt.grid()


plt.subplot(4,1,3)
#plt.plot(motor.t, motor.v_dend, label="Dendrite")
plt.plot(motor.t, motor.v_soma, label="Soma")
#plt.plot(motor.t, motor.v_axon, label="Axon terminal")
#plt.legend()
plt.xlabel("Time (ms)")
plt.ylabel("Membrane Potential (mV)")
plt.title("Motor Neuron")
plt.grid()


plt.subplot(4,1,4)
plt.plot(sensory.t, stim_amp, label="Injected Current")
#plt.legend()
plt.xlabel("Time (ms)")
plt.ylabel("Injected Current (nA)")
plt.title("Injected Stimulus")
plt.grid()
#plt.show()



print("Number of spikes: ", len(list(motor.spike_times)))
