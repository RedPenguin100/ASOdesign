from enum import Enum

import math

import numpy as np
import pandas as pd

from Bio.SeqUtils import gc_fraction

from rna_access import RNAAccess


class GCFilterType(Enum):
    ALL = 'gc'
    TOEHOLD_SITE = 'gc_toehold'
    STEM_SITE = 'gc_stem'


# noinspection DuplicatedCode
class RNATriggerSelector(object):

    # @classmethod
    # def get_triggers(cls, file_name):
    #     df = pd.read_csv(file_name, header=0, index_col=0)
    #     return df

    @classmethod
    def calc_trigger_gc_info(cls, trigger_mrna, trigger_segment_size, trigger_binding_site_size):

        trigger_mrna_size = len(trigger_mrna)
        stem_trigger_site_size = trigger_segment_size - trigger_binding_site_size

        indexes = list(range(0, trigger_mrna_size - trigger_segment_size + 1))
        trigger_segments = [trigger_mrna[i:i + trigger_segment_size] for i in indexes]

        gc_segments = list(map(gc_fraction, trigger_segments))
        gc_segments_bind_site = list(map(lambda x: gc_fraction(x[stem_trigger_site_size:]), trigger_segments))
        gc_segments_stem_site = list(map(lambda x: gc_fraction(x[:stem_trigger_site_size]), trigger_segments))

        d = {
            'trigger_seq': trigger_segments,
            'gc': gc_segments,
            'gc_toehold': gc_segments_bind_site,
            'gc_stem': gc_segments_stem_site,
        }

        df = pd.DataFrame(d, index=indexes)

        return df

    @classmethod
    def calc_access_energies(
            cls, mrna_seq, segment_size, binding_site_size, seed_size, max_span, uuid_str):

        assert binding_site_size >= seed_size

        trigger_mrna_size = len(mrna_seq)
        stem_trigger_site_size = segment_size - binding_site_size

        super_seed_sizes = [seed_size * m for m in range(1, 4)]
        super_seed_sizes = list(filter(lambda x: x <= binding_site_size, super_seed_sizes))
        assert super_seed_sizes
        cls.rna_access = RNAAccess(super_seed_sizes, max_span)
        ra = cls.rna_access
        ra.set_uuid_for_web(uuid_str)
        access_query = [('trigger', mrna_seq)]
        res = ra.calculate(access_query)
        access_res = res['trigger']

        ind_info_list = []
        for pos in range(0, trigger_mrna_size - segment_size + 1):
            pos_info = {}
            for super_seed_size in super_seed_sizes:
                step = super_seed_size // 2
                # n_samples = trigger_binding_site_size / step
                rel_offsets = list(range(0, binding_site_size - super_seed_size, step))
                rel_offsets.append(binding_site_size - super_seed_size)

                abs_offsets = list(map(lambda x: x + pos, rel_offsets))
                bind_energies = access_res[super_seed_size][abs_offsets]
                norm_factor = binding_site_size / super_seed_size
                norm_bind_energies = bind_energies * norm_factor

                # fix last weight relatively to the overlap with the one before
                n_values = math.ceil((binding_site_size - super_seed_size) / step) + 1
                assert len(rel_offsets) == n_values
                weights = [1.0] * n_values
                if len(weights) > 1:
                    last_weight = (rel_offsets[-1] - rel_offsets[-2]) / step
                    weights[-1] = last_weight

                fixed_weight_energies = (np.array(norm_bind_energies) * np.array(weights))
                avg_energy = np.sum(fixed_weight_energies) / np.sum(weights)
                min_energy = np.min(fixed_weight_energies)
                avg_id = f"{super_seed_size}_avg"
                min_id = f"{super_seed_size}_min"
                pos_info.update({avg_id: avg_energy, min_id: min_energy})

            ind_info_list.append((pos, pos_info))

        indexes = list(zip(*ind_info_list))[0]
        records = list(zip(*ind_info_list))[1]
        df = pd.DataFrame(records, index=indexes)

        return df

    @classmethod
    def extract_triggers(
            cls, trigger_mrna, n_trigger_segments, trigger_size, trigger_binding_site_size,
            min_gc, max_gc, gc_ranges, gc_filter_type, temperature,
            access_win_size=120, access_seed_size=6, min_peak_dist=5,
            access_cutoff_ratio=0.25, bind_cutoff_ratio=0.66, w_access=1.0, w_bind=1.0, is_cached=False,
            uuid_str=None):
        """
        :param trigger_mrna: trigger sequence
        :param n_trigger_segments: if positive int it is the number of the best candidates to return for each gc range
        :param trigger_size: trigger size
        :param trigger_binding_site_size: the toehold binding site excluding the stem size binding site
        :param min_gc: min gc ratio integer between 0 and 100
        :param max_gc: max gc ratio integer between 0 and 100
        :param gc_ranges: if not 1 it should be integer of gc sub ranges between min_gc to max_gc
        :param gc_filter_type: GCFilterType Enum to filter by all trigger, the toehold site only or stem site
        :param temperature: temperature for bind energies calculations
        :param access_win_size: the sliding window r_access will use for seeking folding interactions
        :param access_seed_size: the seed size we simulate to check for accessibility segments of multiples of this size
        :param min_peak_dist: returned candidates should have offset non overlapping with this distance
        :param access_cutoff_ratio: if 1.0 does not filter otherwise if between 0 and 1 will filter that access cutoff
        :param bind_cutoff_ratio: if 1.0 does not filter otherwise if between 0 and 1 will filter that bind cutoff
        :param w_access: when generating global score this is its relative weight value
        :param w_bind: when generating global score this is its relative weight value
        :param is_cached: if True try to use cache files to load pre-calculated tables
        :param uuid_str: RNAAccess module create temporal file so need uuid prefix for parallel run
        """
        assert len(trigger_mrna) > 1
        assert len(trigger_mrna) >= trigger_size

        assert trigger_binding_site_size < trigger_size

        assert access_seed_size <= trigger_binding_site_size

        trigger_mrna = trigger_mrna.upper().replace('T', 'U')

        # gc filter
        gc_filter_col = gc_filter_type.value
        gc_info = cls.calc_trigger_gc_info(trigger_mrna, trigger_size, trigger_binding_site_size)
        gc_info['gc_filter'] = gc_info[gc_filter_col]

        access_energies = cls.calc_access_energies(
            trigger_mrna, trigger_size, trigger_binding_site_size, access_seed_size, access_win_size, uuid_str)

        # ae_col1 = f"{access_seed_size}_avg"
        # ae_col2 = f"{access_seed_size * 2}_avg"
        selected_cols = [f"{access_seed_size * i}_avg" for i in [1, 2]]
        filtered_access_energies = access_energies.loc[:, access_energies.columns.isin(selected_cols)]

        access_energies['avg_access'] = filtered_access_energies.mean(axis=1)

        assert gc_ranges >= 1
        gc_values = np.linspace(min_gc, max_gc, num=gc_ranges + 1)
        gc_values /= 100

        df_list = []
        for i in range(gc_ranges):
            inc = 'left' if (i + 1 < gc_ranges) else 'both'

            min_range_gc = gc_values[i]
            max_range_gc = gc_values[i + 1]
            gc_indexes = gc_info[gc_info['gc_filter'].between(min_range_gc, max_range_gc, inclusive=inc)].index

            gc_access_energies = access_energies.loc[gc_indexes]

            gc_range = f"gc_range_{min_range_gc}_{max_range_gc}"

            gc_access_energies['gc_range'] = gc_range

            df_list.append(gc_access_energies)

        df = pd.concat(df_list, join='outer', axis=0).fillna(float('nan'))

        return df


if __name__ == '__main__':
    # from design.toehold.sequence_consts import YEAST_M_CHERRY
    # g_seq = YEAST_M_CHERRY
    # g_n_triggers = -1
    # g_trigger_size = 36
    # g_trigger_binding_site_size = 24
    # g_min_gc = 30
    # g_max_gc = 80
    # g_gc_ranges = 5
    # g_gc_filter = GCFilterType.STEM_SITE
    # g_temperature = 37
    # g_access_win_size = 120
    # g_access_seed_size = 6
    # g_min_peak_dist = 5
    # g_access_cutoff_ratio = 0.25
    # g_bind_cutoff_ratio = 0.66
    # g_w_access = 1.0
    # g_w_bind = 1.0
    # g_is_cached = False

    g_seq = 'AGCCGCUUU'
    g_n_triggers = -1
    g_trigger_size = 6
    g_trigger_binding_site_size = 3
    g_min_gc = 0
    g_max_gc = 100
    g_gc_ranges = 1
    g_gc_filter = GCFilterType.STEM_SITE
    g_temperature = 37
    g_access_win_size = 120
    g_access_seed_size = 3
    g_min_peak_dist = 5
    g_access_cutoff_ratio = 0.25
    g_bind_cutoff_ratio = 0.66
    g_w_access = 1.0
    g_w_bind = 1.0
    g_is_cached = False
    df = RNATriggerSelector.extract_triggers(
        g_seq, g_n_triggers, g_trigger_size, g_trigger_binding_site_size, g_min_gc, g_max_gc, g_gc_ranges,
        g_gc_filter, g_temperature,
        g_access_win_size, g_access_seed_size, g_min_peak_dist,
        g_access_cutoff_ratio, g_bind_cutoff_ratio,
        g_w_access, g_w_bind,
        g_is_cached, )

    print(df)