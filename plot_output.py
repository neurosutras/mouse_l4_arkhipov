import click
import os, sys
from nested.utils import list_find
from bmtk.analyzer.spike_trains import plot_raster


@click.command()
@click.option("--spikes-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True)
@click.option("--with-histogram", type=bool, default=True)
@click.option("--population", type=str, default='l4')
@click.option("--group-by", type=str, default='model_name')
@click.option("--title", type=str, default='Control')
@click.option("--bionet-config-file-path", type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True)
def main(spikes_file_path, with_histogram, population, group_by, title, bionet_config_file_path):
    """
    Generates spike raster plot from L4 simulation output.
    :param spikes_file_path: str (path)
    :param with_histogram: bool
    :param population: str
    :param group_by: str
    :param bionet_config_file_path: str (path)
    """
    plot_raster(config_file=bionet_config_file_path, population=population, with_histogram=with_histogram,
                times=None, title=title, show=True, group_by=group_by, spikes_file=spikes_file_path)


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(os.path.basename(__file__)) != -1, sys.argv) + 1):],
         standalone_mode=False)