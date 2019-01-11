#!/usr/bin/python2
import os
import sys
import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
from simlib import gillespie, run_simulation, save_timeseries_histogram,\
                   save_datafile, save_last_histogram, calculate_results,\
                   seconds2str, plot_avg_vs_time

def main(args):
    """The main method that runs the simulation"""
    ## ----- Configuration

    experiment = {
        "mean_burst_size": args["mean_burst_size"],
        "burst_size_distribution": args["burst_size_distribution"],
        "initial_time": 0,
        "duration": args["duration"],
        "cell_population": args["cell_population"],
        "initial_population": {
            "dna": 1,
            "mrna": 0,
            "mrna1": 0,
            "mrna2": 0,
            "protein": 0,
        },
        "framestep": (float(args["duration"]) /4 if args["framestep"] == 0 \
                      else args["framestep"]),
        "molecule_to_plot": args["molecule_data"],
        "n_jobs": args["n_jobs"],
        "pacifier_active": args["pacifier"],
        "hard_limit": args["hard_limit"],
    }

    if args["silent"]:
        sys.stdout = open(os.devnull, 'w')

    events_description = {
        "burst_arrival": {"rate": args["burst_arrival_rate"], "elem": "dna"},
        "mrna_sene1": {"rate": args["mrna_dis_rate"]*2, "elem": "mrna1"},
        "mrna_decay": {"rate": args["mrna_dis_rate"]*2, "elem": "mrna2"},
        "protein_prod": {"rate": args["protein_prod_rate"], "elem": "mrna"},
        "protein_decay": {"rate": args["protein_dis_rate"], "elem": "protein"},
    }

    exp_id = 'CP%.6d_mbd-%s_mbs%02.2f_kb%.1f_um%.1f_kp%.1f_up%.2f_T%07d_FS%03.2f_i%03d'\
              % (experiment["cell_population"],
                 experiment["burst_size_distribution"],
                 experiment["mean_burst_size"],
                 args["burst_arrival_rate"],
                 args["mrna_dis_rate"],
                 args["protein_prod_rate"],
                 args["protein_dis_rate"],
                 experiment["duration"],
                 experiment["framestep"],
                 args["i_repetition"])

    experiment["exp_id"] = exp_id

    start_time = time.time()

    if os.path.isfile("./data/%s.npy" % experiment["exp_id"]) and not args["redo_simulation"]:
        print("We already have this data. I will plot it")
        data = np.load("./data/%s.npy" % experiment["exp_id"])
        experiment_data = data[2]
        simulation_runned = False
    else:
        try:
            experiment_data = run_simulation(experiment, events_description)
        except (KeyboardInterrupt, SystemExit):
            return 1
        simulation_runned = True

    end_time = time.time()

    if args["plot_timeseries"]:
        save_timeseries_histogram(events_description, experiment, experiment_data)
    if not args["skip_final_plot"]:
        save_last_histogram(events_description, experiment, experiment_data, show_plot=True)
        plot_avg_vs_time(experiment_data)
        plt.show()

    rslt = calculate_results(experiment_data)
    results_str = "% 7d" % rslt[0]
    for i in range(1, len(rslt)):
        results_str += "\t% 15.6f" % rslt[i]
    results_str += "\t% 15.2f" % experiment["mean_burst_size"]
    results_str += "\t%s" % experiment["burst_size_distribution"]
    results_str += "\t%s" % experiment["exp_id"]
    results_str += "\t%s" % seconds2str(end_time - start_time)
    results_str += "\n"


    if args["results_file"] is not None:
        args["results_file"].write(results_str)

    print "\nDone!\n"


parser = argparse.ArgumentParser(description='Gillespi simulation.')

# Main experiment arguments
parser.add_argument('-n', '--n-cell', dest='cell_population', action='store',
                    required=True, type=int,
                    help='The number of cells in the simulation')

parser.add_argument('-d', '--experiment-duration', dest='duration',
                    action='store', required=False, default=1000, type=int,
                    help='The stop time for the simulation')

parser.add_argument('-kb', '--burst-arrival-rate', dest='burst_arrival_rate',
                    action='store', required=False, default=10.0, type=float,
                    help='The rate for mRNA burst arrivals')

parser.add_argument('-um', '--mrna-dis-rate', dest='mrna_dis_rate',
                    action='store', required=False, default=2.0, type=float,
                    help='The rate for mRNA disintegration')

parser.add_argument('-kp', '--protein-prod-rate', dest='protein_prod_rate',
                    action='store', required=False, default=15.0, type=float,
                    help='The rate for protein production')

parser.add_argument('-up', '--protein-dis-rate', dest='protein_dis_rate',
                    action='store', required=False, default=0.05, type=float,
                    help='The rate for protein disintegration')

parser.add_argument('-r', '--redo-simulation', dest='redo_simulation',
                    action='store_true',
                    help=('Force the execution of the simulation, even if '
                          'data is saved for this parameters'))

# mRNA 
parser.add_argument('-mbs', '--mean-burst-size', dest='mean_burst_size',
                    action='store', required=True, type=float,
                    help='The average size of mRNA bursts')

parser.add_argument('-bsd', '--burst-size-distribution',
                    dest='burst_size_distribution', action='store',
                    choices=['geometric', 'conditional_geometric', 'delta'],
                    required=False, default="conditional_geometric", type=str,
                    help='The type of distribution for the mRNA bursts')

# Plots
parser.add_argument('-nfp', '--skip-final-plot', dest='skip_final_plot',
                    action='store_true',
                    help='Do not generate the plot at the end of simulation')

parser.add_argument('--molecule-data', dest='molecule_data', action='store',
                    choices=['protein', 'mrna'],
                    required=False, default="protein", type=str,
                    help='Data to use for the histogram plots')

parser.add_argument('-pts', '--plot-timeseries', dest='plot_timeseries',
                    action='store_true',
                    help='Save PNG sequence of the simulation\'s timeseries')

parser.add_argument('-tsfs', '--timeseries-framestep', dest='framestep',
                    action='store', required=False, default=0, type=float,
                    help='Time intervals between data frames')

parser.add_argument('-rf', '--results-file', dest='results_file',
                    type=argparse.FileType('a'),
                    default=None)

# Multiprocessing jobs and repetitions
parser.add_argument('-nj', '--n-jobs', dest='n_jobs', action='store',
                    required=False, default=4, type=int,
                    help='The number of CPUs to use')

parser.add_argument('-l', '--hard-time-limit', dest='hard_limit', action='store',
                    required=False, default=86400, type=int,
                    help='Stops the simulation after this many seconds')

parser.add_argument('-i', '--i-repetition', dest='i_repetition', action='store',
                    required=False, default=1, type=int,
                    help='The id of experiment repetition')

parser.add_argument('-q', '--quiet', '--no-pacifier', dest='pacifier',
                    action='store_false',
                    help='Supress pacifier output')

parser.add_argument('-s', '--silent', dest='silent', action='store_true',
                    help='Supress all output')

if __name__ == "__main__":
    main(vars(parser.parse_args()))
