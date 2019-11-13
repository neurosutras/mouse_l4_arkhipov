"""
Simulates a model of L4 of V1 consisting of 45,000 individual cells divided into 7 different cell-type classes.
This script interacts with nested.optimize to find a set of parameters controlling synaptic connection strengths for
eachscalar weight vectorsThe goal of this script is get each cell-type to have a population-averaged firing rate set in the optimization configuration file.
n example network of 450 cell receiving two kinds of external input as defined in the configuration file"""
import json
import shutil
import pandas as pd
import click

from nested.parallel import *
from nested.optimize_utils import *
from neuron import h
from collections import defaultdict

from bmtk.simulator import bionet
from bmtk.simulator.bionet.nrn import synaptic_weight
from bmtk.simulator.bionet.modules import SaveSynapses, SpikesMod
from bmtk.analyzer.spike_trains import spike_statistics
from bmtk.utils.io import ioutils


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

context = Context()


@click.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True, ))
@click.option("--config-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False),
              default='config/optimize_L4_bionet_config.yaml')
@click.option("--export", is_flag=True)
@click.option("--output-dir", type=str, default='data')
@click.option("--export-file-path", type=str, default=None)
@click.option("--label", type=str, default=None)
@click.option("--interactive", is_flag=True)
@click.option("--verbose", type=int, default=2)
@click.option("--plot", is_flag=True)
@click.option("--debug", is_flag=True)
@click.option("--simulate", type=bool, default=True)
@click.pass_context
def main(cli, config_file_path, export, output_dir, export_file_path, label, interactive, verbose, plot, debug,
         simulate):
    """

    :param cli: contains unrecognized args as list of str
    :param config_file_path: str (path)
    :param export: bool
    :param output_dir: str
    :param export_file_path: str
    :param label: str
    :param interactive: bool
    :param verbose: int
    :param plot: bool
    :param debug: bool
    :param simulate: bool
    """
    context.update(locals())
    kwargs = get_unknown_click_arg_dict(cli.args)
    context.disp = verbose > 0

    if 'procs_per_worker' not in kwargs:
        kwargs['procs_per_worker'] = int(MPI.COMM_WORLD.size)

    context.interface = get_parallel_interface(source_file=__file__, source_package=__package__, **kwargs)
    context.interface.start(disp=context.disp)
    context.interface.ensure_controller()
    config_optimize_interactive(__file__, config_file_path=config_file_path, output_dir=output_dir, export=export,
                                export_file_path=export_file_path, label=label, disp=context.disp,
                                interface=context.interface, verbose=verbose, plot=plot, debug=debug, **kwargs)
    if simulate:
        run_tests()

    if context.plot:
        context.interface.apply(plt.show)

    if not interactive:
        context.interface.stop()


def config_worker():

    if 'plot' not in context():
        context.plot = False
    if 'verbose' not in context():
        context.verbose = 1
    else:
        context.verbose = int(context.verbose)
    if 'debug' not in context():
        context.debug = False
    start_time = time.time()

    temp_output_dir = context.temp_output_path.replace('.hdf5', '')
    temp_output_dir = context.comm.bcast(temp_output_dir, root=0)
    if context.comm.rank == 0:
        os.mkdir(temp_output_dir)

    ioutils.set_world_comm(context.comm)
    if context.comm != ioutils.bmtk_world_comm.comm:
        raise RuntimeError('optimize_L4_bionet: problem setting bmtk_world_comm')

    if context.verbose > 1 and context.debug:
        print('optimize_L4_bionet: comm rank/size: %i/%i; bionet.io_tools.io.mpi_size: %i' %
          (context.comm.rank, context.comm.size, bionet.io_tools.io.mpi_size))
    sys.stdout.flush()

    bionet_config_dict = json.load(open(context.bionet_config_file_path, 'r'))
    bionet_config_dict['manifest']['$OUTPUT_DIR'] = temp_output_dir
    conf = bionet.Config.from_dict(bionet_config_dict, validate=True)
    conf.build_env()

    spikes_file_path = conf.output['output_dir'] + '/' + conf.output['spikes_file']

    # load network using config.json params
    graph = bionet.BioNetwork.from_config(conf)

    # run initial simulation using config.json params
    sim = bionet.BioSimulator.from_config(conf, network=graph)

    if context.verbose > 1 and context.comm.rank == 0:
        print('optimize_L4_bionet: initialization took %.2f s' % (time.time() - start_time))

    context.update(locals())


def run_tests():
    features = context.interface.execute(compute_features, context.x0_array, context.export)
    sys.stdout.flush()
    time.sleep(1.)
    if len(features) > 0:
        features, objectives = context.interface.execute(get_objectives, features)
    else:
        objectives = dict()
    sys.stdout.flush()
    time.sleep(1.)
    if context.export:
        collect_and_merge_temp_output(context.interface, context.export_file_path, verbose=context.disp)
    sys.stdout.flush()
    time.sleep(1.)
    print('params:')
    pprint.pprint(context.x0_dict)
    print('features:')
    pprint.pprint(features)
    print('objectives:')
    pprint.pprint(objectives)
    sys.stdout.flush()
    time.sleep(1.)
    context.interface.execute(shutdown_worker)
    context.update(locals())


def update_context(x, local_context=None):
    """

    :param x: array
    :param local_context: :class:'Context'
    """
    if local_context is None:
        local_context = context
    x_dict = param_array_to_dict(x, local_context.param_names)

    # update_syn_weights(graph, gradients)


def compute_features(x, export=False):
    """

    :param x: array
    :param export: bool
    :return: dict
    """
    update_source_contexts(x, context)

    conf = context.conf
    sim_step = bionet.BioSimulator(network=context.graph, dt=conf.dt, tstop=conf.tstop, v_init=conf.v_init,
                                   celsius=conf.celsius, nsteps_block=conf.block_step)

    # Attach mod to simulation that will be used to keep track of spikes.

    spikes_recorder = \
        SpikesMod(spikes_file=context.spikes_file_path, tmp_dir='', spikes_sort_order='gid', mode='w')
    sim_step.add_mod(spikes_recorder)

    # run simulation
    sim_step.run()

    # Get the spike statistics of the output, using "groupby" will get averaged firing rates across each model
    spike_stats_df = spike_statistics(context.spikes_file_path, simulation=sim_step, groupby='model_name',
                                      populations='l4')

    results = dict()
    for pop_name, rate_val in spike_stats_df['firing_rate']['mean'].items():
        feature_name = 'mean_rate_' + pop_name
        results[feature_name] = rate_val

    if export:
        updated_weights_dir = context.export_file_path.replace('.hdf5', '')
        os.mkdir(updated_weights_dir)
        connection_recorder = SaveSynapses(updated_weights_dir)
        connection_recorder.initialize(sim_step)
        connection_recorder.finalize(sim_step)
        context.comm.barrier()

    return results


def get_objectives(features, export=False):
    """

    :param features: dict
    :param export: bool
    :return: tuple of dict
    """
    if context.comm.rank == 0:
        objectives = {}
        for objective_name in context.objective_names:
            objectives[objective_name] = ((context.target_val[objective_name] - features[objective_name]) /
                                          context.target_range[objective_name]) ** 2.
        return features, objectives


def update_syn_weights(net, gradients):
    """Go through each cell and update their synaptic weights by the gradient (determined by the cell's model_name).

    If the incoming synapse is excitatory we update the weight by +gradient and if it is inhibitory we update the weight
    by -gradient.
    """

    for gid, cell in net.get_local_cells().items():
        trg_pop = cell['model_name']

        for con in cell.connections():
            if con.is_virtual:
                continue

            src_node = con.source_node
            src_type = src_node['ei']
            con.syn_weight += gradients[trg_pop]*(-1.0 if src_type == 'i' else 1.0)

            con.syn_weight = max(con.syn_weight, 0.0)


def shutdown_worker():
    if os.path.isdir(context.temp_output_dir):
        shutil.rmtree(context.temp_output_dir)


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(os.path.basename(__file__)) != -1, sys.argv) + 1):],
         standalone_mode=False)
