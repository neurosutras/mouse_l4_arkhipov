from bmtk.analyzer.visualization.spikes import plot_raster
import matplotlib.pyplot as plt
import click


@click.command()
@click.option("--spikes-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True)
@click.option("--with-histogram", type=bool, default=True)
@click.option("--with-labels", multiple=True, type=str, default=['l4'])
@click.option("--group-by", type=str, default='model_name')
@click.option("--nodes-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True)
@click.option("--node-types-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True)
def main(spikes_file_path, with_histogram, with_labels, group_by, nodes_file_path, node_types_file_path):
    """
    Generates spike raster plot from L4 simulation output.
    :param spikes_file_path: str (path)
    :param with_histogram: bool
    :param with_labels: list of str
    :param group_by: str
    :param nodes_file_path: str (path)
    :param node_types_file_path: str (path)
    """
    plt.figure('Raster')
    plot_raster(spikes_file_path, with_histogram=with_histogram, with_labels=with_labels, group_by=group_by,
                show_plot=False, nodes_file=nodes_file_path, node_types_file=node_types_file_path)
    plt.show()
