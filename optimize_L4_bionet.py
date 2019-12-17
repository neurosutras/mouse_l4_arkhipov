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
from distutils.dir_util import copy_tree

from bmtk.simulator import bionet
from bmtk.simulator.bionet.nrn import synaptic_weight
from bmtk.simulator.bionet.modules import SaveSynapses, SpikesMod
from bmtk.utils.io import ioutils
from bmtk.utils.reports import SpikeTrains


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

context = Context()


@click.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True, ))
@click.option("--config-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False),
              default='config/optimize_L4_bionet_toy_config.yaml')
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
        if not os.path.isdir(context.output_dir):
            os.mkdir(context.output_dir)
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
    if 'input_dir' in context() and os.path.isdir(context.input_dir):
        bionet_config_dict['manifest']['$INPUT_DIR'] = context.input_dir
    if 'network_dir' in context() and os.path.isdir(context.network_dir):
        bionet_config_dict['manifest']['$NETWORK_DIR'] = context.network_dir
    conf = bionet.Config.from_dict(bionet_config_dict, validate=True)
    conf.build_env()

    spikes_file_path = conf.output['output_dir'] + '/' + conf.output['spikes_file']

    # construct network using config.json params
    graph = bionet.BioNetwork.from_config(conf)

    # initialize simulation using config.json params
    sim = bionet.BioSimulator.from_config(conf, network=graph)

    init_weights = defaultdict(dict)
    for target_gid, cell in graph.get_local_cells().items():
        target_pop_name = cell['model_name']

        for con in cell.connections():
            init_weights[target_pop_name][con] = con.syn_weight

    if context.debug:
        syn_count = defaultdict(lambda: defaultdict(int))
        for target_pop_name in init_weights:
            for con in init_weights[target_pop_name]:
                if con.is_virtual:
                    source_pop_name = con.source_node._population
                else:
                    source_pop_name = con.source_node['model_name']
                syn_count[target_pop_name][source_pop_name] += 1
        syn_count = defaultdict_to_dict(syn_count)
        syn_count_list = context.comm.gather(syn_count, root=0)
        if context.comm.rank == 0:
            gathered_syn_count = defaultdict(lambda: defaultdict(int))
            for this_syn_count in syn_count_list:
                for target_pop_name in this_syn_count:
                    for source_pop_name in this_syn_count[target_pop_name]:
                        gathered_syn_count[target_pop_name][source_pop_name] += \
                            this_syn_count[target_pop_name][source_pop_name]
            print('synapse counts per projection')
            pprint.pprint(defaultdict_to_dict(gathered_syn_count))
            sys.stdout.flush()
            time.sleep(.1)
        context.comm.barrier()

    if context.comm.rank == 0:
        node_ids = {'grating_rate': [], 'grating_pref_rate': [], 'grating_ortho_rate': []}
        graph_node_props_df = graph.get_node_groups(populations='l4').groupby('model_name')
        for (pop_name, node_props_df) in graph_node_props_df:
            (i, first_row) = next(node_props_df.iterrows())
            if first_row['ei'] == 'i':
                node_ids['grating_rate'].extend(node_props_df.node_id.values)
            else:
                for category, target_angle in [('grating_pref_rate', context.grating_angle),
                                               ('grating_ortho_rate', (context.grating_angle + 90.) % 360.)]:
                    node_props_df_copy = node_props_df.copy()
                    node_props_df_copy['sort_val'] = abs(node_props_df_copy.tuning_angle - target_angle)
                    node_props_df_copy = node_props_df_copy.sort_values('sort_val')
                    node_ids[category].extend(node_props_df_copy.node_id.values[:50])

        epochs = {'spont_rate': {'start': 0.,
                                 'stop': 500.,},
                  'grating_rate': {'start': 500.,
                                   'stop': 1000.,
                                   'node_ids': node_ids['grating_rate']},
                  'grating_pref_rate': {'start': 500.,
                                        'stop': 1000.,
                                        'node_ids': node_ids['grating_pref_rate']},
                  'grating_ortho_rate': {'start': 500.,
                                        'stop': 1000.,
                                        'node_ids': node_ids['grating_ortho_rate']}
                  }

        if context.debug:
            print('Firing rate analysis epochs:')
            pprint.pprint(epochs)
            sys.stdout.flush()
            time.sleep(.1)
    else:
        epochs = None
    epochs = context.comm.bcast(epochs, root=0)

    if context.verbose > 0 and context.comm.rank == 0:
        print('optimize_L4_bionet: initialization took %.2f s' % (time.time() - start_time))
        sys.stdout.flush()
        time.sleep(.1)

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

    local_context.weight_factors = defaultdict(dict)

    for param_name, param_val in x_dict.items():
        param_type, target_pop_name, pop_category = param_name.split('.')
        if param_type == 'weight_factor':
            for source_pop_name in context.pop_categories[pop_category]:
                local_context.weight_factors[target_pop_name][source_pop_name] = param_val

    scale_projection_weights(graph=local_context.graph, weight_factors=local_context.weight_factors,
                             init_weights=local_context.init_weights)

    if local_context.debug and local_context.comm.rank == 0:
        print('weight_factors:')
        pprint.pprint(local_context.weight_factors)
        sys.stdout.flush()
        time.sleep(.1)


def compute_features(x, export=False):
    """

    :param x: array
    :param export: bool
    :return: dict
    """
    update_source_contexts(x, context)
    start_time = time.time()

    conf = context.conf
    sim_step = bionet.BioSimulator(network=context.graph, dt=conf.dt, tstop=conf.tstop, v_init=conf.v_init,
                                   celsius=conf.celsius, nsteps_block=conf.block_step)

    if context.export:
        # Attach mod to simulation that will record and cache spikes to disk and export to sonata format (h5).
        spikes_recorder = \
            SpikesMod(spikes_file=context.conf.output['spikes_file'], tmp_dir=context.temp_output_dir,
                      spikes_sort_order='gid', mode='w')
        sim_step.add_mod(spikes_recorder)
    else:
        # Record spikes but do not cache to disk
        sim_step.set_spikes_recording()

    # run simulation
    sim_step.run()

    if context.verbose > 1 and context.comm.rank == 0:
        print('optimize_L4_bionet: pid: %i; simulation with x: %s took %.2f s' %
              (os.getpid(), str(list(x)), time.time() - start_time))
        sys.stdout.flush()
        time.sleep(.1)

    if export:
        export_dir = context.export_file_path.replace('.hdf5', '')
        if context.comm.rank == 0:
            if not os.path.isdir(export_dir):
                os.mkdir(export_dir)
        context.comm.barrier()

    if export and not context.debug:
        start_time = time.time()
        connection_recorder = SaveSynapses(export_dir)
        connection_recorder.initialize(sim_step)
        connection_recorder.finalize(sim_step)
        if context.verbose > 1 and context.comm.rank == 0:
            print('optimize_L4_bionet: pid: %i; exporting connection weights with x: %s took %.2f s' %
                  (os.getpid(), str(list(x)), time.time() - start_time))
            sys.stdout.flush()
            time.sleep(.1)
        context.comm.barrier()

    start_time = time.time()
    # Get the average firing rates per cell type in each epoch
    if export:
        if context.comm.rank == 0:
            firing_rates_dict = \
                get_firing_rates_by_cell_type_from_file(context.spikes_file_path, simulation=sim_step,
                                                        population='l4', epochs=context.epochs)
    else:
        local_firing_rates_dict = \
            get_local_firing_rates_by_cell_type_from_sim(sim_step, population='l4', epochs=context.epochs)
        list_of_local_firing_rates_dict = context.comm.gather(local_firing_rates_dict, root=0)
        if context.comm.rank == 0:
            firing_rates_dict = dict()
            for local_firing_rates_dict in list_of_local_firing_rates_dict:
                for epoch_name in local_firing_rates_dict:
                    if epoch_name not in firing_rates_dict:
                        firing_rates_dict[epoch_name] = dict()
                    for cell_type in local_firing_rates_dict[epoch_name]:
                        if cell_type not in firing_rates_dict[epoch_name]:
                            firing_rates_dict[epoch_name][cell_type] = []
                        firing_rates_dict[epoch_name][cell_type].extend(local_firing_rates_dict[epoch_name][cell_type])
    context.comm.barrier()

    if context.comm.rank == 0:
        if export:
            copy_tree(context.temp_output_dir, export_dir)

        for epoch_name in firing_rates_dict:
            for cell_type in firing_rates_dict[epoch_name]:
                firing_rates_dict[epoch_name][cell_type] = np.mean(firing_rates_dict[epoch_name][cell_type])

        results = dict()
        for epoch_name in firing_rates_dict:
            for pop_name, rate_val in firing_rates_dict[epoch_name].items():
                feature_name = '%s.%s' % (epoch_name, pop_name)
                results[feature_name] = rate_val

        if context.verbose > 1:
            print('optimize_L4_bionet: pid: %i; analysis with x: %s took %.2f s' %
                  (os.getpid(), str(list(x)), time.time() - start_time))
            sys.stdout.flush()
            time.sleep(.1)

        if context.debug:
            pprint.pprint(firing_rates_dict)
            sys.stdout.flush()
            time.sleep(.1)
            context.update(locals())

        return results


def get_firing_rates_by_cell_type_from_file(spikes_file, simulation, population='l4', cell_type_attr_name='model_name',
                                            epochs=None):
    """

    :param spikes_file: path to .h5 file
    :param simulation: :class:'BioSimulator'
    :param populations: str
    :param cell_type_attr_name: str
    :param epochs: dict: {str: tuple of float (ms)}
    :return: dict
    """
    def get_firing_rates_by_cell_type(spike_trains, population, node_ids, cell_type_dict, start, stop):
        """

        """
        rate_dict = dict()
        duration = (stop - start) / 1000.  # sec
        for gid in node_ids:
            cell_type = cell_type_dict[gid]
            if cell_type not in rate_dict:
                rate_dict[cell_type] = []
            spike_times = spike_trains.get_times(gid, population=population, time_window=[start, stop])
            rate = len(spike_times) / duration  # Hz
            rate_dict[cell_type].append(rate)
        return rate_dict

    spike_trains = SpikeTrains.load(spikes_file)
    cell_type_dict = dict()
    all_nodes_df = simulation.net.get_node_groups(populations=population)
    all_node_ids = all_nodes_df['node_id'].values
    for gid, cell_type in all_nodes_df[['node_id', cell_type_attr_name]].values:
        cell_type_dict[gid] = cell_type
    sim_end = simulation.simulation_time(units='ms')

    rate_dict = dict()
    if epochs is None:
        rate_dict['all'] = \
            get_firing_rates_by_cell_type(spike_trains, population, all_node_ids, cell_type_dict, start=0.,
                                          stop=sim_end)
    else:
        for epoch_name, epoch_dict in epochs.items():
            if 'node_ids' in epoch_dict and len(epoch_dict['node_ids']) > 0:
                this_node_ids = epoch_dict['node_ids']
            else:
                this_node_ids = all_node_ids
            rate_dict[epoch_name] = \
                get_firing_rates_by_cell_type(spike_trains, population, this_node_ids, cell_type_dict,
                                              start=epoch_dict['start'], stop=epoch_dict['stop'])

    return rate_dict


def get_local_firing_rates_by_cell_type_from_sim(simulation, population='l4', cell_type_attr_name='model_name',
                                                 epochs=None):
    """

    :param simulation: :class:'BioSimulator'
    :param populations: str
    :param cell_type_attr_name: str
    :param epochs: dict: {str: tuple of float (ms)}
    :return: dict
    """
    def get_firing_rates_by_cell_type(spikes_table, node_ids, cell_type_dict, start, stop):
        """

        """
        rate_dict = dict()
        duration = (stop - start) / 1000.  # sec
        for gid in node_ids:
            cell_type = cell_type_dict[gid]
            if cell_type not in rate_dict:
                rate_dict[cell_type] = []
            spike_times = np.array(spikes_table[gid])
            spike_times = spike_times[(start <= spike_times) & (spike_times <= stop)]
            rate = len(spike_times) / duration  # Hz
            rate_dict[cell_type].append(rate)
        return rate_dict

    spikes_table = simulation.spikes_table
    local_cells = simulation.net.get_local_cells()
    gid_pool = simulation.net.gid_pool
    cell_type_dict = dict()
    all_local_node_ids = []

    for gid, cell in local_cells.items():
        if population == gid_pool.get_pool_id(gid).population:
            cell_type = cell[cell_type_attr_name]
            cell_type_dict[gid] = cell_type
            all_local_node_ids.append(gid)
    sim_end = simulation.simulation_time(units='ms')

    rate_dict = dict()
    if epochs is None:
        rate_dict['all'] = \
            get_firing_rates_by_cell_type(spikes_table, all_local_node_ids, cell_type_dict, start=0., stop=sim_end)
    else:
        for epoch_name, epoch_dict in epochs.items():
            if 'node_ids' in epoch_dict and len(epoch_dict['node_ids']) > 0:
                this_node_ids = set(epoch_dict['node_ids']).intersection(all_local_node_ids)
            else:
                this_node_ids = all_local_node_ids
            rate_dict[epoch_name] = \
                get_firing_rates_by_cell_type(spikes_table, this_node_ids, cell_type_dict, start=epoch_dict['start'],
                                              stop=epoch_dict['stop'])

    return rate_dict


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


def scale_projection_weights(graph, weight_factors, init_weights=None):
    """
    Iterate over local cells and connections, and scale weights according to target and source population. If
    init_weights are specified, the initial weights for each connection are scaled. Otherwise, the current weight value
    of each connection is scaled.
    :param graph: :class:'BioNetwork'
    :param weight_factors: dict: {target_pop_name (str): {source_pop_name (str): float} }
    :param init_weights: dict: {target_gid (int): {connection (:class:'ConnectionStruct'): float} }
    """
    for target_gid, cell in graph.get_local_cells().items():
        target_pop_name = cell['model_name']

        for con in cell.connections():
            if con.is_virtual:
                source_pop_name = con.source_node._population
            else:
                source_pop_name = con.source_node['model_name']

            if target_pop_name not in weight_factors or source_pop_name not in weight_factors[target_pop_name]:
                continue

            if init_weights is not None:
                if con not in init_weights[target_pop_name]:
                    raise KeyError('scale_projection_weights: no initial weight stored for connection between %s cell '
                                   '%i and %s cell %i' %
                                   (target_pop_name, target_gid, source_pop_name, con.source_node.gid))
                init_weight = init_weights[target_pop_name][con]
            else:
                init_weight = con.syn_weight
            con.syn_weight = init_weight * weight_factors[target_pop_name][source_pop_name]


def shutdown_worker():
    if context.comm.rank == 0:
        if os.path.isdir(context.temp_output_dir):
            shutil.rmtree(context.temp_output_dir)


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(os.path.basename(__file__)) != -1, sys.argv) + 1):],
         standalone_mode=False)
