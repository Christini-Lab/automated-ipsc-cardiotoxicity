import myokit
import matplotlib.pyplot as plt
import numpy as np


proto = myokit.Protocol()

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))

for v in range(-90, 60, 10):
    proto.add_step(-100, 1000)
    proto.add_step(v, 1000)

for rseries in [2E-3, 3E-3, 10E-3, 20E-3]:
    mod = myokit.load_model('./mmt-files/kernik_2019_NaL_art.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    rs = mod.get('artifact.r_series')
    rs.set_rhs(rseries)

    v_cmd = mod.get('artifact.v_cmd')
    v_cmd.set_rhs(0)
    v_cmd.set_binding('pace') # Bind to the pacing mechanism
    holding_proto = myokit.Protocol()
    holding_proto.add_step(-80, 30000)
    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())
    mem = mod.get('artifact')

    v_cmd = mod.get('artifact.v_cmd')
    v_cmd.set_binding(None)
    v_cmd.set_binding('pace')

    t_max = proto.characteristic_time()
    times = np.arange(0, t_max, 0.1)

    sim = myokit.Simulation(mod, proto)

    dat = sim.run(t_max, log_times=times)

    i_out = [val/60 for val in dat['artifact.i_out']] 
    axs[1].plot(times, i_out)
    axs[2].plot(times, dat['ina.i_Na'], label=rseries)


axs[0].plot(times, dat['artifact.v_cmd'])
axs[0].set_ylabel('Voltage')
axs[1].set_ylabel('I_out')
axs[2].set_ylabel('I_Na')
axs[2].set_xlabel('Time')
axs[2].legend()

plt.show()



