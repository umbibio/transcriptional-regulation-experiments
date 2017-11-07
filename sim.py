#!/usr/bin/python
import argparse, os
import numpy as np
from simlib import gillespie, run_simulation, save_timeseries_histogram,\
                   save_datafile, save_last_histogram, calculate_results

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
            "protein": 0,
        },
        "framestep": args["framestep"],
        "molecule_to_plot": args["molecule_data"],
    }

    events_description = {
        "burst_arrival": {"rate": 1, "elem": "dna"},
        "mrna_decay": {"rate": 2, "elem": "mrna"},
        "protein_prod": {"rate": 15, "elem": "mrna"},
        "protein_decay": {"rate": 0.01, "elem": "protein"},
    }

    exp_id = 'CP%.6d_mbd-%s_mbs%02.2f_kb%.1f_um%.1f_kp%.1f_up%.2f_T%07d_FS%03d'\
              % (experiment["cell_population"],
                 experiment["burst_size_distribution"],
                 experiment["mean_burst_size"],
                 events_description["burst_arrival"]["rate"],
                 events_description["mrna_decay"]["rate"],
                 events_description["protein_prod"]["rate"],
                 events_description["protein_decay"]["rate"],
                 experiment["duration"],
                 experiment["framestep"])
    
    experiment["exp_id"] = exp_id

    if os.path.isfile("../data/%s.npy" % experiment["exp_id"]) and not args["redo_simulation"]:
        print("We already have this data. I will plot it")
        data = np.load("../data/%s.npy" % experiment["exp_id"])
        experiment_data = data[2]
    else:
        experiment_data = run_simulation(experiment, events_description)

    if args["plot_timeseries"]:
        save_timeseries_histogram(events_description, experiment, experiment_data)
    if not args["skip_final_plot"]:
        save_last_histogram(events_description, experiment, experiment_data, show_plot=True)

    rslt = calculate_results(experiment_data)
    results_str = "% 7d" % rslt[0]
    for i in range(1, len(rslt)):
        results_str += "\t% 15.4f" % rslt[i]
    results_str += "\t% 15.4f" % experiment["mean_burst_size"]
    results_str += "\t%s" % experiment["burst_size_distribution"]
    results_str += "\t%s" % experiment["exp_id"]
    results_str += "\n"
    print results_str


    if args["results_file"] is not None:
        args["results_file"].write(results_str)

    print "\nDone!\n"


parser = argparse.ArgumentParser(description='Gillespi simulation.')

# Main experiment arguments
parser.add_argument('--n-cell', dest='cell_population', action='store',
                    required=True, type=int,
                    help='The number of cells in the simulation')

parser.add_argument('--experiment-duration', dest='duration', action='store',
                    required=False, default=1000, type=int,
                    help='The stop time for the simulation')

parser.add_argument('--redo-simulation', dest='redo_simulation', action='store_true',
                    help='Force the execution of the simulation, even if data is saved for this parameters')

# mRNA 
parser.add_argument('--mean-burst-size', dest='mean_burst_size', action='store',
                    required=True, type=float,
                    help='The average size of mRNA bursts')

parser.add_argument('--burst-size-distribution', dest='burst_size_distribution', action='store',
                    choices=['conditional_geometric', 'delta'],
                    required=False, default="conditional_geometric", type=str,
                    help='The type of distribution for the mRNA bursts')

# Plots
parser.add_argument('--skip-final-plot', dest='skip_final_plot', action='store_true',
                    help='Do not generate the plot at the end of simulation')

parser.add_argument('--molecule-data', dest='molecule_data', action='store',
                    choices=['protein', 'mrna'],
                    required=False, default="protein", type=str,
                    help='Data to use for the histogram plots')

parser.add_argument('--plot-timeseries', dest='plot_timeseries', action='store_true',
                    help='Save PNG sequence of the simulation\'s timeseries')

parser.add_argument('--timeseries-framestep', dest='framestep', action='store',
                    required=False, default=50, type=int,
                    help='Time intervals between data frames')

parser.add_argument('--results-file', dest='results_file', type=argparse.FileType('a'),
                    default=None)

main(vars(parser.parse_args()))
