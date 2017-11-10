"""Methods for the simulation"""
import gc, os
from time import time
from multiprocessing import Process, Pool
import numpy as np
import matplotlib.pyplot as plt

def secondsToStr(t):
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll,b : divmod(ll[0],b) + ll[1:],
            [(t*1000,),1000,60,60])

def run_batch(experiment, events_description, njobs, job):
    timeseries = []
    dataseries = []
    pacifier = 0.0
    samples = experiment["cell_population"] / njobs
    pacifier_inc = float(max(100, samples)) / 100

    prefix = "\033[F" * (njobs - job)
    suffix = "\n" * (njobs - job - 1)

    seed = int((time() - int(time()/10**8) * 10**8)) * (job + 1)
    np.random.seed(seed)

    for iteration in range(samples):
        if iteration >= np.ceil(pacifier):
            pacifier += pacifier_inc
            print prefix + "Job #: %s\tProgress: % 4.0f%%" % (job, (float(pacifier)/samples * 100 )) + suffix

        time_list, data_list = gillespie(experiment, events_description)
        timeseries.append(time_list)
        dataseries.append(data_list)
    
    return timeseries, dataseries

def run_simulation(experiment, events_description):

    njobs = 4

    print "\nStarting simulation...\n"
    print "Total cells:\t\t%s" % experiment["cell_population"]
    print "Mean mRNA Burst Size:\t%s" % experiment["mean_burst_size"]
    print "mRNA BurstArrvl rate:\t%s" % (events_description["burst_arrival"]["rate"])
    print "mRNA Dis rate:\t\t%s" % events_description["mrna_decay"]["rate"]
    print "Protein Prd rate:\t%s" % events_description["protein_prod"]["rate"]
    print "Protein Dis rate:\t%s" % events_description["protein_decay"]["rate"]
    print "Time to simulate:\t%s" % experiment["duration"]
    print "\n" * njobs

    # Each element of this list represents a number of proteins created by one mRNA
    #timeseries = []
    #dataseries = []

    pool = Pool(njobs)
    data = {}
    time_s = {}
    data_s = {}

    for job in range(njobs):
        data[job] = pool.apply_async(run_batch, (experiment, events_description, njobs, job))

    pool.close()
    pool.join()

    timeseries, dataseries = [], []

    for job in range(njobs):
        time_s[job], data_s[job] = data[job].get()
        timeseries += time_s[job]
        dataseries += data_s[job]

    print "\nRearranging arrays...\n\n"
    protein_number = []
    mrna_number = []
    time = []
    pacifier = 0.0
    for i in range(int(np.round(experiment["duration"] / experiment["framestep"])) + 1):
        if i >= int(np.ceil(pacifier)):
            pacifier += experiment["duration"] / 100
            print "\033[FProgress: % 4d%%" % int(np.ceil(float(i) / experiment["duration"] * 100 ))

        time.append(timeseries[0][i])
        protein_number.append([])
        mrna_number.append([])
        for cid in range(len(timeseries)):
            protein_number[i].append(dataseries[cid][i]["protein"])
            mrna_number[i].append(dataseries[cid][i]["mrna"])
    print "\033[FProgress: % 4d%%" % 100

    experiment_data = {"time": time, "protein": protein_number, "mrna": mrna_number}

    print "\nSaving datafiles..."
    save_datafile(events_description, experiment, experiment_data)
    return experiment_data

def gillespie(experiment, events_description):

    def calc_effective_rate():
        effective_rate = 0.0
        for event in events_description:
            effective_rate += molecule_population[event_generator[event]] * event_rate[event]
        return effective_rate

    prob_success = {}
    prob_success["conditional_geometric"] = 1.0 / experiment["mean_burst_size"]
    prob_success["geometric"] = 1.0 / (experiment["mean_burst_size"] + 1)

    def geometric_distribution():
        return int(np.random.geometric(prob_success["geometric"]) - 1)

    def conditional_geometric_distribution():
        return np.random.geometric(prob_success["conditional_geometric"])

    def delta_distribution():
        return experiment["mean_burst_size"]

    distributions = {
        "delta": delta_distribution,
        "conditional_geometric": conditional_geometric_distribution,
        "geometric": geometric_distribution,
    }

    clock_time = experiment["initial_time"]
    molecule_population = {}
    event_generator = {}
    event_rate = {}
    event_probability = {}
    event_threshold = {}

    for event_name in events_description:
        elem = events_description[event_name]["elem"]
        molecule_population[elem] = experiment["initial_population"][elem]
        event_generator[event_name] = elem

        rate = float(events_description[event_name]["rate"])
        event_rate[event_name] = rate
        event_probability[event_name] = 0.0
        event_threshold[event_name] = 0.0

    effective_rate = calc_effective_rate()
    timeframe = 0.0
    time = []
    data = []

    events_order = [
        "burst_arrival",
        "mrna_decay",
        "protein_prod",
        "protein_decay"
        ]

    stop_time = experiment["initial_time"] + experiment["duration"]
    while clock_time < stop_time and effective_rate > 0.0:

        while clock_time >= timeframe:
            time.append(timeframe)
            data.append(molecule_population.copy())
            timeframe += experiment["framestep"]

        for event in events_description:
            event_probability[event] = \
                molecule_population[event_generator[event]]\
                * event_rate[event] / effective_rate

        threshold = 0.0
        for event in events_order:
            threshold += event_probability[event]
            event_threshold[event] = threshold

        choose_event = np.random.uniform()
        if choose_event < event_threshold["burst_arrival"]:
            molecule_population["mrna"] += distributions[experiment["burst_size_distribution"]]()
        elif choose_event < event_threshold["mrna_decay"]:
            molecule_population["mrna"] -= 1
        elif choose_event < event_threshold["protein_prod"]:
            molecule_population["protein"] += 1
        elif choose_event < event_threshold["protein_decay"]:
            molecule_population["protein"] -= 1

        clock_time += np.random.exponential(1.0 / effective_rate)
        effective_rate = calc_effective_rate()

    while clock_time >= timeframe:
        time.append(timeframe)
        data.append(molecule_population.copy())
        timeframe += experiment["framestep"]

    time.append(timeframe)
    data.append(molecule_population)
    return time, data

def save_datafile(events_description, experiment, experiment_data):
    filepath = '../data/%s' % (experiment["exp_id"])
    np.save(filepath, [events_description, experiment, experiment_data])

def save_last_histogram(events_description, experiment, experiment_data, max_bin=None, max_height=None, show_plot=False):
    
    final_time = experiment_data["time"][-1]
    final_molecule_n = experiment_data[experiment["molecule_to_plot"]][-1]
    
    mean = np.mean(final_molecule_n)
    max_bin = np.ceil(mean * 4)
    bin_size = max([int(np.ceil(float(max_bin) / 50)), 1])

    number, bins, patches = plt.hist(final_molecule_n, bins=xrange(0, bin_size*50, bin_size), normed=True)
    max_height = np.amax(number) * 1.2
    plt.clf()

    save_histogram(events_description, experiment, final_time, final_molecule_n, max_bin=max_bin, max_height=max_height, show_plot=show_plot)

def save_timeseries_histogram(events_description, experiment, experiment_data, max_bin=None, max_height=None):

    final_molecule_n = experiment_data[experiment["molecule_to_plot"]][-1]
    
    mean = np.mean(final_molecule_n)
    max_bin = np.ceil(mean * 4)
    bin_size = max([int(float(max_bin) / 50), 1])
    number, bins, patches = plt.hist(final_molecule_n, bins=xrange(0, bin_size*50, bin_size), normed=True)
    max_height = np.amax(number) * 1.2
    plt.clf()

    for i, time in enumerate(experiment_data["time"]):
        save_histogram(events_description, experiment, time, experiment_data[experiment["molecule_to_plot"]][i], max_bin=max_bin, max_height=max_height, make_own_dir=True)
        plt.clf()


def save_histogram(events_description, experiment, clock_time, protein_numbers, max_bin=None, max_height=None, make_own_dir=False, show_plot=False):

    p1_mean = np.mean(protein_numbers)
    p1_var = np.var(protein_numbers)
    if p1_mean > 0.0:
        p1_ff = p1_var / p1_mean
    else:
        p1_ff = 0.0

    if max_bin is None:
        bin_size = int(np.ceil(p1_mean / 20))
    else:
        bin_size = max([int(float(max_bin) / 50), 1])

    fig1 = plt.figure(1)

    number, bins, patches = plt.hist(protein_numbers,
                                     bins=xrange(0, bin_size*50, bin_size),
                                     normed=True,
                                     label='Simulated distribution')

    ax = plt.gca()
    if max_height is not None:
        ax.set_ylim([0.0, max_height])
    else:
        ax.set_ylim([0.0, np.amax(number) * 1.2])
    plt.xlabel('%s number' % experiment["molecule_to_plot"])
    plt.ylabel('Probability')
    title = ax.set_title("%d Cells %s distribution\ndistribution of $m_b$: %s\n"
              "$<m_b>=%s$, $k_b=%.2f$, $\mu_m=%s$, $k_p=%s$, $\mu_p=%s$\n"
              "$t=%.1f$, $<N>=%.2f$, $\sigma^2=%.2f$, $FF=%.2f$"
              % (experiment["cell_population"],
                 experiment["molecule_to_plot"],
                 experiment["burst_size_distribution"],
                 experiment["mean_burst_size"],
                 events_description["burst_arrival"]["rate"],
                 events_description["mrna_decay"]["rate"],
                 events_description["protein_prod"]["rate"],
                 events_description["protein_decay"]["rate"],
                 clock_time,
                 p1_mean, p1_var, p1_ff))
    title.set_y(1.05)
    fig1.subplots_adjust(top=0.75)
    plt.grid(True)
    plt.legend()

    if make_own_dir:
        directory = '../figs/tseries/%s/%s' % (experiment["exp_id"], experiment["molecule_to_plot"])
        if not os.path.exists(directory):
            os.makedirs(directory)
        fig_path = directory + '/T%07d.png' % (clock_time)
    else:
        fig_path = '../figs/%s_%s_T%07d.png' % (experiment["molecule_to_plot"], experiment["exp_id"], clock_time)

    print "\rSaving fig to: %s" % fig_path
    fig1.savefig(fig_path)

    gc.collect()

def plot_avg_vs_time(experiment_data):
    avgs = {}
    for key in experiment_data:
        if key <> "time":
            avgs[key] = []

    for i, time in enumerate(experiment_data["time"]):
        for molecule in avgs:
            avgs[molecule].append(np.mean(experiment_data[molecule][i]))

    # Two subplots, the axes array is 1-d
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].set_title('Avg Populations vs time')
    axarr[-1].set_xlabel("Time")
    for i, key in enumerate(avgs):
        axarr[i].plot(experiment_data["time"], avgs[key])
        axarr[i].set_ylabel("Avg %s #" % key)

def calculate_results(experiment_data):
    final_protein_numbers = experiment_data["protein"][-1]
    final_mrna_numbers = experiment_data["mrna"][-1]

    total_cells = len(final_mrna_numbers)

    protein_mean = np.mean(final_protein_numbers)
    protein_var = np.var(final_protein_numbers)
    protein_ff = protein_var / protein_mean
    mrna_mean = np.mean(final_mrna_numbers)
    mrna_var = np.var(final_mrna_numbers)
    mrna_ff = mrna_var / mrna_mean

    return [total_cells, mrna_mean, mrna_var, mrna_ff, protein_mean, protein_var, protein_ff]