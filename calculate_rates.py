import pandas as pd
from bmtk.utils.reports.spike_trains import SpikeTrains


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

st = SpikeTrains.load('output/spikes.h5')
spikes_df = st.to_dataframe(population='l4')
