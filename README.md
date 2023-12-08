# RaVeN
Non-parametric inference of the radial velocity distribution of stars in globular clusters.
RaVeN relies on [FIGARO](https://github.com/sterinaldi/FIGARO) for the non-parametric inference, so please be sure that you have it installed and properly working (pip installation is sometimes faulty...).

If you use RaVeN in your research, please cite [Ramírez-Tannus & Rinaldi (2024)](https://uncyclopedia.com/wiki/Frankly_Disappointing_Telescope).

© 2023 Stefano Rinaldi, María Claudia Ramírez-Tannus

## Basic usage
Install RaVeN with `python setup.py install` (add `--user` if you don't have admin privileges on your machine). The folder structure that RaVeN expects is the following:

```
cluster
└── stars
    ├─ star_1.txt
    ├─ star_2.txt
    └─ star_3.txt
```

The naming convention does not matter as long as the folder structure is respected. Each .txt file must contain two columns with radial velocity measurement and error for each available epoch (0 if no error is available).

To analyse a cluster with RaVeN, simply run the command line instruction `raven -i path/to/cluster/stars`. If you want to specify some bounds for the output plot, add the option `-b "[vmin, vmax]"` (remember the quotation marks). A different number of draws can be set with `-n` and the non-parametric inference, once performed at least one time, can be skipped with `-p`. The available options can be displayed with `raven -h`. After the run, the folder will look like this:

```
cluster
├── stars
│   ├─ star_1.txt
│   ├─ star_2.txt
│   └─ star_3.txt
├── draws
│   ├─ star_1.pkl
│   ├─ star_2.pkl
│   ├─ star_3.pkl
│   ├─ draws_cluster.pkl
│   └─ posteriors_single_event.pkl
├── events
│   ├─ star_1.txt
│   ├─ star_2.txt
│   └─ star_3.txt
├── cluster.pdf
├── prob_cluster.txt
├── no_errors.txt
├── single_fraction_cluster.pdf
└── p_single_cluster.txt
```

Among these files, `cluster.pdf` reports the radial velocity probability density, `single_fraction_cluster.pdf` contains the posterior probability density for the fraction of single stars and `p_single_cluster.txt` contains the objects included in the cluster ranked by their probability of being single stars. If some of the measurements do not have an associated uncertainty, `no_errors.txt` will be created stating the relative uncertainty that has been associated with the measurements.
