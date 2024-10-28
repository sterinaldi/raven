# RAVEN – RAdial VElocity, Non-parametric
Non-parametric inference of the radial velocity distribution of stars and identification of binaries in clusters. 

If you use RAVEN in your research, please cite [Rinaldi & Ramírez-Tannus (2024)](https://uncyclopedia.com/wiki/Frankly_Disappointing_Telescope). The main result from this paper can be reproduced using the content of this repository: if you are interested in doing so, please find the instructions at the bottom of this page.

## Basic usage
Install RAVEN with 
```
git clone git@github.com:sterinaldi/raven.git
cd raven
pip install .
```

The folder structure expected by RAVEN is the following:

```
cluster
└── stars
    ├─ star_1.txt
    ├─ star_2.txt
    └─ star_3.txt
```

The naming convention does not matter as long as the folder structure is respected. Each `.txt` file must contain two columns with radial velocity measurement and error for each available epoch (0 if no error is available).

To analyse a cluster with RAVEN, simply run the command line instruction `raven -i path/to/cluster/stars`. If you want to specify some bounds for the output plot, add the option `-b "[vmin, vmax]"` (remember the quotation marks). A different number of draws can be set with `-n` and the non-parametric inference, once performed at least one time, can be skipped with `-p`. The available options can be displayed with `raven -h`. After the run, the folder will look like this:

```
cluster
├── stars
│   ├─ star_1.txt
│   ├─ star_2.txt
│   └─ star_3.txt
├── draws
│   ├─ star_1.json
│   ├─ star_2.json
│   ├─ star_3.json
│   ├─ draws_cluster.json
│   └─ posteriors_single_event.json
├── events
│   ├─ star_1.txt
│   ├─ star_2.txt
│   └─ star_3.txt
├── cluster.pdf
├── prob_cluster.txt
├── no_errors.txt
├── single_fraction_cluster.pdf
├── samples_fraction_cluster.txt
├── p_single_cluster.pdf
└── p_single_cluster.txt
```

### Output files overview
Content of the output files produced by RAVEN:

* `cluster.pdf`: radial velocity probability density;
* `single_fraction_cluster.pdf`: posterior probability density for the fraction of single stars;
* `samples_fraction_cluster.txt`: posterior samples for the fraction of single stars;
* `p_single_cluster.txt`: objects included in the cluster ranked by their probability of being single stars;
* `p_single_cluster.pdf`: fancy plot displaying the probability of each object of being a single star;
* `no_errors.txt`: if some of the measurements do not have an associated uncertainty, this file contains the relative uncertainty that has been associated with the measurements.

## M17
[Rinaldi & Ramírez-Tannus (2024)](https://uncyclopedia.com/wiki/Frankly_Disappointing_Telescope) presents the analysis of 20 O-type stars in the giant star-forming region M17 reported in  [Ramírez-Tannus et al. (2024)](https://www.aanda.org/articles/aa/full_html/2024/10/aa50256-24/aa50256-24.html). Our results can be reproduced using the code and data stored in this repository.

Once RAVEN is installed, move to the `M17` directory and run the analysis:
```
cd M17
raven -i stars -b "[-40,60]" --sana_variability
```
The analysis should finish in around 10 minutes on a normal laptop. The plots included in the paper are `M17.pdf` and `p_single_M17.pdf`.
