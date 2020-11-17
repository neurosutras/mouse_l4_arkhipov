import json
import click
import os, sys
import matplotlib.pyplot as plt
from nested.utils import list_find
from bmtk.simulator import bionet
from bmtk.analyzer.spike_trains import plot_raster


@click.command()
@click.option("--spikes-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True)
@click.option("--with-histogram", type=bool, default=True)
@click.option("--population", type=str, default='l4')
@click.option("--group-by", type=str, default='model_name')
@click.option("--title", type=str, default='Control')
@click.option("--bionet-config-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True)
@click.option("--nodes-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True)
@click.option("--node-types-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True)
def main(spikes_file_path, with_histogram, population, group_by, title, bionet_config_file_path, nodes_file_path,
         node_types_file_path):
    """
    Generates spike raster plot from L4 simulation output.
    :param spikes_file_path: str (path)
    :param with_histogram: bool
    :param population: str
    :param group_by: str
    :param bionet_config_file_path: str (path)
    :param nodes_file_path: str (path)
    :param node_types_file_path: str (path)
    """
    plot_raster(config_file=bionet_config_file_path, population=population, with_histogram=with_histogram,
                times=None, title=title, show=True, group_by=group_by, spikes_file=spikes_file_path,
                nodes_file=nodes_file_path, node_types_file=node_types_file_path)
    # plot_raster(spikes_file_path, with_histogram=with_histogram, with_labels=with_labels, group_by=group_by,
    #            show_plot=False, nodes_file=nodes_file_path, node_types_file=node_types_file_path)


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(os.path.basename(__file__)) != -1, sys.argv) + 1):],
         standalone_mode=False)