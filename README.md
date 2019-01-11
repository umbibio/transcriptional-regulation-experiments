```
$ ./sim.py --help
usage: sim.py [-h] -n CELL_POPULATION [-d DURATION] [-kb BURST_ARRIVAL_RATE]
              [-um MRNA_DIS_RATE] [-kp PROTEIN_PROD_RATE]
              [-up PROTEIN_DIS_RATE] [-r] -mbs MEAN_BURST_SIZE
              [-bsd {geometric,conditional_geometric,delta}] [-nfp]
              [--molecule-data {protein,mrna}] [-pts] [-tsfs FRAMESTEP]
              [-rf RESULTS_FILE] [-nj N_JOBS] [-l HARD_LIMIT]
              [-i I_REPETITION] [-q] [-s]

Gillespi simulation.

optional arguments:
  -h, --help            show this help message and exit
  -n CELL_POPULATION, --n-cell CELL_POPULATION
                        The number of cells in the simulation
  -d DURATION, --experiment-duration DURATION
                        The stop time for the simulation
  -kb BURST_ARRIVAL_RATE, --burst-arrival-rate BURST_ARRIVAL_RATE
                        The rate for mRNA burst arrivals
  -um MRNA_DIS_RATE, --mrna-dis-rate MRNA_DIS_RATE
                        The rate for mRNA disintegration
  -kp PROTEIN_PROD_RATE, --protein-prod-rate PROTEIN_PROD_RATE
                        The rate for protein production
  -up PROTEIN_DIS_RATE, --protein-dis-rate PROTEIN_DIS_RATE
                        The rate for protein disintegration
  -r, --redo-simulation
                        Force the execution of the simulation, even if data is
                        saved for this parameters
  -mbs MEAN_BURST_SIZE, --mean-burst-size MEAN_BURST_SIZE
                        The average size of mRNA bursts
  -bsd {geometric,conditional_geometric,delta}, --burst-size-distribution {geometric,conditional_geometric,delta}
                        The type of distribution for the mRNA bursts
  -nfp, --skip-final-plot
                        Do not generate the plot at the end of simulation
  --molecule-data {protein,mrna}
                        Data to use for the histogram plots
  -pts, --plot-timeseries
                        Save PNG sequence of the simulation's timeseries
  -tsfs FRAMESTEP, --timeseries-framestep FRAMESTEP
                        Time intervals between data frames
  -rf RESULTS_FILE, --results-file RESULTS_FILE
  -nj N_JOBS, --n-jobs N_JOBS
                        The number of CPUs to use
  -l HARD_LIMIT, --hard-time-limit HARD_LIMIT
                        Stops the simulation after this many seconds
  -i I_REPETITION, --i-repetition I_REPETITION
                        The id of experiment repetition
  -q, --quiet, --no-progress
                        Supress progress output
  -s, --silent          Supress all output
```