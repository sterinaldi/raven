# RaVeN
Non-parametric inference of the radial velocity distribution of stars in globular clusters.
RaVeN relies on [FIGARO](https://github.com/sterinaldi/FIGARO) for the non-parametric inference, so please be sure that you have it installed and properly working (pip installation is sometimes faulty...).

Basic usage: install RaVeN with `python setup.py install` (add `--user` if you don't have admin privileges on your machine). The folder structure that RaVeN expects is the following:

```
cluster
└── stars
    ├─ star_1.txt
    ├─ star_2.txt
    └─ star_3.txt
```

The naming convention does not matter as long as the folder structure is respected. Each .txt file must contain two columns with radial velocity measurement and error for each available epoch.

To analyse a cluster with RaVeN, simply run the command line instruction `raven -i path/to/cluster/stars`. If you want to specify some bounds for the output plot, add the option `-b "[vmin, vmax]"` (remember the quotation marks). A different number of draws can be set with `-n` and the non-parametric inference, once performed once, can be skipped with `-p`. The available options can be displayed with `raven -h`. After the run, the folder will look like this:

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
├── log_cluster.pdf
├── prob_cluster.txt
└── probability_single_star.txt
```

Among the other files, the most useful ones are `cluster.pdf`, that reports the radial velocity probability density, and `probability_single_star.txt`, which contains the objects included in the cluster sorted by their probability of being single stars.

If you use RaVeN in your research, please cite [Ramírez-Tannus & Rinaldi (2024)](https://uncyclopedia.com/wiki/Frankly_Disappointing_Telescope).

© 2023 Ste Rinaldi, María Claudia Ramírez-Tannus
