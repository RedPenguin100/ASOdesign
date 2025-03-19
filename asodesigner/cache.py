import os
import pickle

from asodesigner.consts import CACHE_DIR


def load_cache_off_target_hybridization(organism):
    CACHE_DIR.mkdir(exist_ok=True)
    cache = f'off-target-cache-hybridization-{organism}.pickle'
    cache_path = CACHE_DIR / cache

    if os.path.exists(str(cache_path)):
        with open(str(cache_path), 'rb') as file:
            loaded_data = pickle.load(file)
    else:
        loaded_data = dict()

    return loaded_data, cache_path


def load_cache_off_target_wc(organism):
    CACHE_DIR.mkdir(exist_ok=True)
    cache = f'off-target-cache-wc-{organism}.pickle'
    cache_path = CACHE_DIR / cache

    if os.path.exists(cache_path):
        with open(cache_path, 'rb') as file:
            loaded_data = pickle.load(file)
    else:
        loaded_data = dict()

    return loaded_data, cache_path


def save_cache(cache_path, loaded_data):
    with open(str(cache_path), 'wb') as file:
        pickle.dump(loaded_data, file)


def update_loaded_data(loaded_data, results_dict):
    for result in results_dict:
        if result not in loaded_data:
            loaded_data[result] = results_dict[result]
