# Hierarchical Object Search on Spatial Networks

This repository contains a framework to evaluate the performance of the COL-Tree index for object search on spatial networks. The data structures were first presented in a research paper at ICAPS and subsequently extended for a journal paper (under review). The repository provides the experimental framework, implementation of competing algorithms, and scripts to reproduce all results in the papers. We have released it under the permissive MIT license. Please cite the ICAPS paper if you use our code:

Tenindra Abeywickrama, Muhammad Aamir Cheema, and Sabine Storandt. [Hierarchical Graph Traversal for Aggregate k Nearest Neighbors Search in Road Networks](https://doi.org/10.1609/icaps.v30i1.6639). In ICAPS 2020, 2-10.

```
@inproceedings{abey2020hierarchical,
  title     = {Hierarchical Graph Traversal for Aggregate k Nearest Neighbors Search in Road Networks},
  author    = {Abeywickrama, Tenindra and Cheema, Muhammad Aamir and Storandt, Sabine},
  booktitle = {Proceedings of the 30th International Conference on Automated Planning and Scheduling},
  year      = {2020},
  pages     = {2--10},
  doi       = {10.1609/icaps.v30i1.6639},
  url       = {https://doi.org/10.1609/icaps.v30i1.6639}
}
```

# Requirements

## Dependencies

Before the executables can be compiled, the following packages/libraries must be installed.

- g++ with C++23 support (tested on version 15.2.0)
- Boost.Serialization 1.70 or newer (tested on version 1.90)
- CMake 3.20 or newer (3.21 or newer to use the build presets below)
- gnuplot for figure generation
- [Google Sparsehash](http://code.google.com/p/sparsehash/)
- [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/download)
- Optional: OpenMP for parallel processing of untimed index construction

Note: On Ubuntu, we installed the libraries with `apt install libboost-serialization-dev libsparsehash-dev libmetis-dev`. CMake checks for all three libraries when configuring and reports any that are missing.

## System & Hardware

We ran experiments on a system with Ubuntu 26.04, an 8-core/16-thread CPU, 64GB RAM, and a 1TB SSD. A system with less RAM (e.g., 32GB) should be able to run USA-t experiments and less storage may be sufficient if not running exhaustive parameter testing experiments.

# Directory Structure

The C++ source code can be found in the top-level directory of the repository, and in the subdirectories `command`, `processing`, `queue`, `tuple`, and `utility`. All bash scripts for various tasks, such as setting up the directories, running experiments, and generating figures, are found in the `lmk_scripts` subdirectory.


# Compilation

1. From the top-level directory of the repository (the one containing `CMakeLists.txt`), generate makefiles in a `build` directory using `CMake`. The build defaults to a `Release` (`-O3 -march=native`) build optimized for the CPU of the build machine. To use a specific `g++` version, add `-DCMAKE_CXX_COMPILER=/path/to/g++` to this command.

    ```
    cmake -S . -B build
    ```

2. Compile all executables.

    ```
    cmake --build build -j
    ```

    This produces three executables in `build`: `rn_obj_search` (used for query times), `rn_obj_search_stats` (collects additional query statistics) and `rn_obj_search_par` (enables parallelization for certain untimed tasks).

Alternatively, with CMake 3.21 or newer, use the presets defined in `CMakePresets.json`. The `release` preset is equivalent to the commands above, and the `debug` preset builds with debug symbols and compiler warnings into `build-debug`.

```
cmake --preset release
cmake --build --preset release -j
```

Note: The build directory records absolute paths, so if the repository directory is moved or renamed, delete the build directory and repeat the steps above.


# Setup

1. Open the `globalVariables` bash script in the `lmk_scripts` directory in your favorite editor.

2. Change the `output_path` variable to the full path of an **existing** (i.e., create it first) location where you wish to store all data (e.g., indexes, object sets, figures).

    Note: This path **must** be different to the path containing the code.

3. Change the `exe_path` variable to the full path of the build directory created above (e.g., `/home/user/rn-obj-search/build`, if the repository was cloned to `/home/user/rn-obj-search`).

4. Run `resetExperimentalSetup` in the `lmk_scripts` directory to create all necessary sub-directories.

    ```
    cd lmk_scripts
    bash resetExperimentalSetup
    ```

5. While in the `lmk_scripts` directory, run `downloadZenodo` to download the datasets from [Zenodo](https://doi.org/10.5281/zenodo.23533198) to the `output_path` directory. Zenodo stores the dataset archive `spatial_network_for_object_search.tar.gz` in parts. The script downloads the parts, checks each against its MD5 checksum on Zenodo, joins them into the archive, unpacks it in `output_path` and checks the unpacked files against the archive's `SHA256SUMS` file. It needs `curl` and `python3`. If it stops part-way, run it again: it keeps the parts already downloaded.

    ```
    bash downloadZenodo
    ```

    Note: To download the parts by hand instead, save every file named `spatial_network_for_object_search.tar.gz.<number>.part` from Zenodo to the `output_path` directory and run the commands below in that directory.

    ```
    cat spatial_network_for_object_search.tar.gz.*.part > spatial_network_for_object_search.tar.gz
    tar -zxvf spatial_network_for_object_search.tar.gz
    sha256sum -c SHA256SUMS
    ```

6. Check that `output_path` now contains two directories: `output_path/data` and `output_path/real_world_pois`. The former contains the graph and coordinates files of the networks listed under Datasets below. `output_path/real_world_pois` contains a subdirectory with real-world POI sets for the NW and USA road networks. The archive `spatial_network_for_object_search.tar.gz` (and its parts, if downloaded by hand) is no longer needed and can be deleted.

7. **Optional:** If you would like to run experiments for other [DIMACS](http://www.dis.uniroma1.it/challenge9/) road network datasets, create a new directory `output_path/data/dimacs`, and download the relevant DIMACS graph files and coordinate files for each road network. Add a dataset config for each network (see Datasets), then use the command below to preprocess all DIMACS datasets in this directory. Note that you will need to either disable real-world POI experiments (experiments on synthetic POI sets will still work) or create corresponding real-world POI datasets yourself. If you would like to use completely different road networks (e.g., from OpenStreetMap) or other spatial networks, these need to be processed similarly to ensure the pipeline works (check `transformInputData`).

    ```
    bash transformInputData
    ```

# Datasets

Each dataset has a config file in `lmk_scripts/datasets` named after it with its edge weight as a suffix (`-t` travel time, `-d` distance), for example the main dataset USA-t. A config defines the dataset's index, query, and experimental parameters, and other applicable properties. To add a dataset, copy a config and edit it. Before running experiments, binary graph files must be created for all datasets. First change directory into the `lmk_scripts` directory and then run:

```
bash buildBinaryGraphs
```

Experiments in the following section are run per graph dataset, which is set by changing `default_network` in `lmk_scripts/globalVariables` or by passing `-d <dataset>` to scripts that accept it. As the entire workflow for running the USA dataset used in the paper is quite long, we included the smaller NW dataset so that it can be used to quickly test whether the pipeline is functioning. The table below lists all the datasets included in the Zenodo repository with matching config files.

| Name | Vertices | Edges | Edge-Types | Notes |
|------|----------|-------|------------|-------|
| USA | 23,947,347 | 57,708,624 | t | Road network (DIMACS USA). All query types |
| W | 6,262,104 | 15,119,284 | t, d | Road network (DIMACS W, western USA). AkNN |
| NW | 1,089,933 | 2,545,844 | t, d | Road network (DIMACS NW), small, for quick tests. `t`: all query types; `d`: AkNN |
| NACA0015 | 1,039,183 | 6,229,636 | d | Mesh: the DIMACS10 NACA0015 airfoil mesh. FkN |
| N39W107.SRTMGL3 | 1,442,401 | 11,524,800 | d | Terrain: the SRTM 3-arcsecond elevation tile N39W107 (1,201 × 1,201 cells). AkNN |


# Running Experiments

The following instructions describe the step-by-step procedure to reproduce all figures for the AkNN query results from the paper (followed by how to adapt the procedure for other query experiments and other datasets). First, ensure you are in the `lmk_scripts` directory:

1. Generate sets of query locations for all experiments. The query and object sets of steps 1 and 2 are drawn from `set_seed` in `lmk_scripts/globalVariables`, so the same seed generates the same sets.

    ```
    bash generateQuerySets
    ```

2. Generate the synthetic object (POI) sets. Note: Running this script again for the same dataset with another `set_seed` invalidates any object index built in steps 3 and 4, so the script refuses to replace existing sets unless `-o` is given (`bash generateObjectSets -o`). If you replace the existing object sets, follow steps 3 and 4 again with `-i objects`, which rebuilds only the object indexes.

    ```
    bash generateObjectSets
    ```

3. Build the indexes that the dataset's query types share: the road network indexes (PHL and ALT) and the object indexes built from them and the object sets (R-trees and NVDs), for both synthetic and real-world POI sets (the latter only if applicable). `dataset_indexes` in `globalVariables` sets which of them are built, and the dataset's config sets their parameters (`dataset_index_parameters`). This step is needed once per dataset; `-i` rebuilds only the indexes it lists (e.g., `bash buildDatasetIndexes -i nvd`), and `-i objects` rebuilds the enabled object indexes without the road network indexes.

    ```
    bash buildDatasetIndexes
    ```

    Note: PHL construction is slow on the terrain, mesh and distance-weighted datasets. See Q&A 2 for a faster parallel build (`phl=p`).

4. Build the indexes of the query type given by `-e`, here the SUL-Tree and COL-Trees for the AkNN experiment (add `-o` to rebuild them, and/or `-i objects` to only build the COL-Trees).

    ```
    bash buildExperimentIndexes -e aknn
    ```

5. Run AkNN query time and query statistics experiments and generate figures.

    ```
    bash runQueryExperiments time -e aknn
    bash runQueryExperiments stats -e aknn
    ```

    Note: The above commands also create figures and output them to the `output_path/fig` directory. The raw data is stored in the `output_path/stats` directory. You can recreate figures without running experiments again with `-f`, e.g., `bash runQueryExperiments time -e aknn -f`.

    Note: The query experiment scripts (`run*QueryExperiments` and the query steps of `runParameterExperiments`) pin the process to the CPU core given by `pin_core` in `globalVariables` using `taskset` (Linux, util-linux), so that every run times its queries on the same core. Unpinned, the scheduler picks the core for each run, and cores that service frequent device interrupts (e.g., a GPU's display interrupts) run queries slower. Choose a core without busy interrupts (see `/proc/interrupts`) whose SMT sibling, if any, is idle (see `lscpu -e`), or leave `pin_core` empty to run unpinned. Do not run two timed scripts at the same time, as they will share the same core due to the pinning.

6. To run experiments for other queries, repeat step 4 with the query type given with `-e` ("knn", "aknn", "fkn" or "range"), then run its experiments as in step 5. Step 4 builds the indexes for the parameters of that query type, which sometimes differ (an existing partition is used, so only a query type with another fanout or maximum leaf size, e.g., FkN on USA, builds a new one). The query types share the indexes of step 3, so it is not repeated. Do not repeat step 2 for the same dataset: regenerating the object sets with another `set_seed` invalidates every object index.

    ```
    bash buildExperimentIndexes -e fkn
    bash runQueryExperiments time -e fkn
    bash runQueryExperiments stats -e fkn
    ```

7. The following sequence of commands can be used to reproduce the parameter testing figures in the paper using the `runParameterExperiments` script. Re-run the last two commands to reproduce the results for queries other than AkNN by changing the mode (e.g., `build_fkn_obj` and `fkn`). The first command (to build the network distance indexes) only needs to be run once. The parameter testing SUL-Trees and COL-Trees are written to `output_path/parameter_testing/indexes` and `output_path/parameter_testing/obj_indexes` (`par_dir` in `globalVariables`), separate from the experiment indexes; the PHL and ALT indexes of step 3 and the object sets of step 2 are shared with the experiments.

    ```
    bash runParameterExperiments build_rn
    bash runParameterExperiments build_aknn_obj
    bash runParameterExperiments aknn
    ```

    `bash runParameterExperiments build_rn` uses `rn_obj_search_par` to parallelize SUL-Tree construction using OpenMP. Available threads are split between building multiple SUL-Trees and processing each SUL-Tree's nodes based on available system memory (the split is printed). Set `OMP_NUM_THREADS` to use fewer cores (e.g., `OMP_NUM_THREADS=4 bash runParameterExperiments build_rn`). Refer to the code for further details. Note that all timings reported in the paper are single-threaded; this parallelization is simply for convenience as parameter testing can take quite a long time.

8. `bash runSDLExperiment` compares three ways of building subgraph distance lists (SDLs): Dijkstra with a hash table, Dijkstra without one, and border-vertex-restricted Dijkstra without one. For each method, it builds the SUL-Tree in memory with the dataset's SDL parameters (`sdl_index_parameters` in its config). It then reports the SDL computation time, search spillover γ and index size, plus the index size it would have without the subgraph vertex ordering. Results go to `output_path/stats/<network>_sdl_stats.txt`. It uses the step 4 partition with the same fanout and leaf size, building it first if needed.

9. To reproduce figures for the other datasets included in the paper or for the quick test via the smaller NW dataset, change `default_network` in `lmk_scripts/globalVariables` or pass `-d <dataset>` to scripts that accept it. First execute steps 1 and 2, which are required for each dataset (see note below re: datasets that have different edge-weights like W-t and W-d). Then do step 3, and steps 4 and 5 for each query type in the dataset's config (e.g., N39W107.SRTMGL3 runs AkNN only), as in step 6.

    Note: For datasets that differ only in edge weights (i.e., they have the same graph structure), such as W-t and W-d, the paper's results run the queries on the same query and object set vertices. To do the same, instead of generating new query and object sets in steps 1 and 2, pass `-c <dataset>` to `generateQuerySets` and `generateObjectSets` (e.g., `bash generateQuerySets -d W-d -c W-t` and `bash generateObjectSets -d W-d -c W-t`). The copy fails if the two graphs do not have the same vertices and neighbors.


# Verify Methods

The `runVerification` script checks the results of the query methods against ground truth. It covers one query type, given with `-e` ("knn", "aknn", "fkn" or "range"), and the synthetic object sets. The ground truth is computed with one full Dijkstra search per query vertex (per query vertex set for AkNN) and saved, so the methods can be re-verified without recomputing it. Verification keeps its own copies of the object sets and object indexes in `output_path/verification/obj_indexes` and the ground truth in `output_path/verification/ground_truth`, and does not touch the experiment indexes or statistics. Steps 1 to 4 of Running Experiments must have been completed for that query type (the verification copies the object sets of step 2 and the object indexes of steps 3 and 4). From the `lmk_scripts` directory:

1. Run all verification steps for a query type, e.g., AkNN (whose indexes were built in steps 3 and 4).

    ```
    bash runVerification all -e aknn
    ```

    This runs `objects` (copies the synthetic object sets and their object indexes of the types to `verify_object_indexes` in `globalVariables`, built by `buildDatasetIndexes` and `buildExperimentIndexes`, which must have been already run for the query type), `ground_truth` (computes and saves the ground truth for the query sets of the experiment, in parallel with OpenMP if available; `verify_threads` limits the number of threads) and `verify` (runs every method enabled for the experiment in `globalVariables`, with the queries in parallel for all methods but the COL-Tree and G-tree ones, whose indexes keep per-query state, and reports the number of failing queries per method; if any query fails, it also lists any object indexes of the failing methods that have different objects to the matching object set (to sanity check mismatching), and the script exits with an error).

2. The steps can also be run individually, for example to re-verify after a code change without recomputing the ground truth, or to verify specific methods only with `-m` (`-m all` verifies every method enabled in `globalVariables` either through experiment default methods, e.g., `aknn_methods`, or `verify_methods` if set).

    ```
    bash runVerification verify -e aknn
    bash runVerification verify -e aknn -m colt_phl,ier_phl
    ```

3. To print a saved ground truth file (the queries and their results):

    ```
    bash runVerification dump -e aknn -g <output_path>/verification/ground_truth/NW-t_aknn_varyingk_max.gt
    ```

4. To verify other query types, run `bash runVerification all -e <query type>` for each of them (after building its indexes, see step 6 above). The verification object indexes and ground truth files of different query types coexist, as their file names carry the index parameters and the query type. After `bash generateObjectSets` regenerates the object sets and steps 3 and 4 rebuild their object indexes, `bash runVerification all -e <query type>` must be re-run as well, so that the verification object sets, object indexes and ground truth come from the same object sets.

Note: On the NW dataset the whole verification takes seconds. On the USA dataset the ground truth step is the slow one, as it runs one Dijkstra search over the entire road network per distinct query vertex, and the verify step of the slowest methods (IER-PHL and NVD-PHL for AkNN queries) takes a while even with its queries in parallel. `-d` selects a dataset other than `default_network`, e.g., `bash runVerification all -e aknn -d USA-t`.


# Q&A

1. **Is there a debug mode or can I enable asserts?**

    The `rn_obj_search` and `rn_obj_search_par` executables include some basic error checking that will exit after printing an error. Additionally, the Debug build (`cmake --preset debug`) and `rn_obj_search_stats` include asserts that run further checks and print `Debug check ... failed`, after their file, line and function, and then abort.


2. **PHL construction is extremely slow on distance edge-weights. What can I do?**

    We optimized the original PHL implementation by replacing its priority queue with a custom radix heap and by using AVX-512 instructions (when the CPU supports them) for some subroutines. Both are on by default and cut construction time by about 10%. For more speed, set `phl=p` (instead of `phl=1`) in `dataset_indexes` in `globalVariables`. This builds PHL with `rn_obj_search_par`, using one thread per core unless `OMP_NUM_THREADS` sets the number. We saw roughly a 4.5x speed-up for NACA0015 on our 8C/16T CPU. The labels are still exact but slightly larger (2-4%), which affects distance computation for all methods equally, and the construction time is printed rather than recorded in the index stats. The parallel build also needs more memory (see Q&A 3). It is only a convenience and does not affect the paper's conclusions. Alternatively, consider shortest path indexes that outperform PHL on distance edge-weights, such as [H2CL](https://github.com/henningkoehlernz/road-networks).
    

3. **The parallel PHL build (`phl=p`) runs out of memory. What can I do?**

    The parallel build needs up to about 1.7 times the memory of the sequential one, and PHL labels are large on the terrain, mesh and distance-weighted datasets. For example, the sequential build for N39W107.SRTMGL3-d peaked at about 38GB (about 2.4 hours on our machine), while the parallel build exceeded our 64GB. You can build sequentially (`phl=1`), use fewer threads (e.g., `OMP_NUM_THREADS=4 bash buildDatasetIndexes -i phl`), or limit malloc arenas (`GLIBC_TUNABLES=glibc.malloc.arena_max=2 bash buildDatasetIndexes`). The last two options cost some time.

4. **I don't have enough memory to run the experiments for NACA0015 and/or N39W107.SRTMGL3. Are there any other options?**

    The `transformInputData` command converts the entirety of the SRTM elevation tiles in `output_path/data/srtm` (e.g., `N39W107.SRTMGL3.hgt.zip`) and the DIMACS10 meshes in `output_path/data/dimacs10` (e.g., `NACA0015.graph.bz2`, with the coordinates in `num-coords.tar.bz2`) into graphs by default. However, it is also possible to convert only a crop using the `input_crops` variable in `globalVariables`. Some sample crops are already added in `input_crops` (although the SRTM crop is for a more granular input file, N39W107.SRTMGL1, from the 1-arcsecond SRTM dataset).

5. **How could SMT and PBO affect the experiments?**

    The query experiments time each query on a single thread, pinned to the core given by `pin_core` (see step 5 of Running Experiments). With simultaneous multithreading (SMT, or Hyper-Threading on Intel CPUs) enabled, work could be allocated to the pinned core's sibling (listed by `lscpu -e`) and slow queries down as they share the same physical core. Therefore, it may be best to disable SMT before running any timed experiments (either in BIOS or via software if your operating system supports it). AMD's Precision Boost Overdrive (PBO) raises the CPU's power limits and can raise its maximum boost clock beyond its specification. Since this is unpredictable by nature, it may be best to disable PBO in BIOS.

6. **Why does `localsearch-3` (or `tracker-miner-fs-3`) use CPU and disk while the experiments run?**

    GNOME's file indexer, LocalSearch (formerly Tracker), indexes the files in the folders it is set to search, which can be the whole home directory (as on our machine). If `output_path` is in such a folder, it also indexes the index, object set and statistics files that the experiments write, taking CPU time and disk access away from them. To avoid any impact on the experiments, it may be best to create a file named `.trackerignore` in `output_path` or in a directory above it, as the indexer skips any directory that contains such a file, e.g., if `output_path` is in `~/experiments`:

    ```
    touch ~/experiments/.trackerignore
    ```