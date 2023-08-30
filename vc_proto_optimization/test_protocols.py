import myokit
import matplotlib.pyplot as plt
import numpy as np


def plot_all_currs():
    fig, axs = plt.subplots(5, 1, sharex=True, figsize=(12, 8))

    proto = myokit.Protocol()
    #proto.add_step(-80, 500)

    for v in [-20, .1, 20]:
        proto.add_step(-80, 1000)
        proto.add_step(-100, 30)
        proto.add_step(v, 50)

    for v in [-10, 5, 20]:
        proto.add_step(-80, 500)
        proto.add_step(-45, 400)
        proto.add_step(v, 100)

    for v in [-50, -35, -20]:
        proto.add_step(-80, 500)
        proto.add_step(30, 300)
        proto.add_step(v, 10)
        proto.add_step(10, 50)

    mod = myokit.load_model('./mmt-files/kernik_2019_NaL_art.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    rs = mod.get('artifact.r_series')
    rseries = 5E-3
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
    axs[0].plot(times, dat['artifact.v_cmd'])
    axs[1].plot(times, i_out)
    axs[2].plot(times, dat['ina.i_Na'], label=rseries)
    axs[3].plot(times, dat['ikr.i_Kr'], label=rseries)
    axs[4].plot(times, dat['ical.i_CaL'], label=rseries)


    axs[0].set_ylabel('Voltage')
    axs[1].set_ylabel('I_out')
    axs[2].set_ylabel('I_Na')
    axs[3].set_ylabel('I_Kr')
    axs[4].set_ylabel('I_CaL')
    axs[4].set_xlabel('Time')

    plt.show()


def plot_ikr():
    fig, axs = plt.subplots(4, 1, sharex=True, figsize=(12, 8))

    proto = myokit.Protocol()

    for v in [-40, -40, -40]:
        proto.add_step(-80, 500)
        proto.add_step(20, 300)
        proto.add_step(v, 10)
        proto.add_step(10, 50)

    for v in [-40, -40, -40]:
        proto.add_step(-80, 500)
        proto.add_step(20, 300)
        proto.add_step(v, 10)
        proto.add_step(10, 50)

    mod = myokit.load_model('./mmt-files/kernik_2019_NaL_art.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    rs = mod.get('artifact.r_series')
    rseries = 5E-3
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
    axs[0].plot(times, dat['artifact.v_cmd'])
    axs[1].plot(times, i_out)
    axs[2].plot(times, dat['ina.i_Na'], label=rseries)
    axs[3].plot(times, dat['ikr.i_Kr'], label=rseries)


    axs[0].set_ylabel('Voltage')
    axs[1].set_ylabel('I_out')
    axs[2].set_ylabel('I_Na')
    axs[3].set_ylabel('I_Kr')
    axs[3].set_xlabel('Time')

    plt.show()

plot_all_currs()
