# Mouse V1 Layer 4 network

## Installation
Install NEURON (7.4+) with the python API.

Get the latest version of bmtk from github and install it:
```bash
$ git clone https://github.com/AllenInstitute/bmtk.git
$ cd bmtk
$ pip install -e .  # or python setup.py develop
```

## Running the (toy) network.
We've included a toy example of the L4 network - uses the same cell/synaptic models but only significantly reduced in size to run and debug on a laptop.
To setup the simulation:
 1. First need to use NEURON to compile special mechanisms used by our Allen Cell Typed Database models:
 ```bash
 $ cd biophys_components/mechanisms
 $ nrnivmodl modfiles/
 ```

 2. Run the simulation, either on a single core:
 ```bash
 $ python run run_bionet.py config.json
 ```
 or using multiple cores:
 ```bash
 $ mpirun -n {N} nrniv -mpi -python run_bionet.py
 ```

The simulator will load up the network, then run a few iterations each time adjusting the synaptic weights. The main output is the network
spikes of the last simulation which is stored in **output/spikes.h5**. You can use the script **plot_output.py** to show the results or use
h5py/HDFView also.


## Running the Full L4 network
First [download the layer4 network files from dropbox](https://www.dropbox.com/sh/rfgqv9oqgu4vei9/AAA_nbTeKjU0UjTNmt5alS3za/network?dl=0&subfolder_nav_tracking=1)
 into the **network/** directory. Then open up **config.json** and in the _manifest_ section change ```$NETWORK_DIR``` value from ```$BASE_DIR/toy_network``` to ```$BASE_DIR/network```.

Next go to the [inputs directory](https://www.dropbox.com/sh/rfgqv9oqgu4vei9/AAD3ldnZJMRS23I4ObCyk9sPa/inputs?dl=0&subfolder_nav_tracking=1)
 and copy the files into **inputs/**. In the config under the _inputs_ section change the ```input_file``` values for both the lgn and tw (traveling wave) inputs

```$INPUT_DIR/lgn_spikes.h5``` --> ```$INPUT_DIR/LGN_spikes.g217.t0.h5```

and

```$INPUT_DIR/tw_spikes.h5``` --> ```$INPUT_DIR/TW_spikes.g217.t0.h5```

Running the simulation should be the same. The time it takes to complete will depend on the cluster you're using but for a 1 second simulation we recommend at minimum 20 cores with 64GB.


## Modifying the simulation(s)

**config.json** contains parameters for loading and running the (initial) simulation.

**run_bionet.py** contains the main code for loading the network, running multiple simulations in serial, each time updating the synaptic weights. During the last simulation the network
output (spikes) and updated synaptic weights are saved in the **output** directory.
