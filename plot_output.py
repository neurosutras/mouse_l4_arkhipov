from bmtk.analyzer.visualization.spikes import plot_raster
import matplotlib.pyplot as plt

# Raster plot of the v1 spikes.
plt.figure('Raster')
plot_raster('output/spikes.h5', with_histogram=True, with_labels=['l4'], group_by='model_name', show_plot=False,
            nodes_file='toy_network/l4_nodes.h5', node_types_file='toy_network/l4_node_types.csv')
plt.show()
