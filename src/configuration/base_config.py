import copy

BASE_CONFIG = {
    "experiment_name": "base",
    "data": {},
    "model": {"layer_sizes": [29260, 50, 10, 1], "dropout_p": 0.3},
    "epochs": 25,
    "lr": 0.0001,
    "loss": "mse",
    "checkpoint_dir": "checkpoints",

}

def apply_diff(base_cfg, diff):
    cfg = copy.deepcopy(base_cfg)
    for k, v in diff.items():
        if isinstance(v, dict):
            cfg[k] = apply_diff(base_cfg[k], v)
        else:
            cfg[k] = v
    return cfg

def first_run():
    diffs = [
        {"experiment_name": "base"},
        {"epochs": 50,
         "data": {"batch_size":16},
         "experiment_name": "bs_16"},
        {"epochs": 50,
         "data": {"batch_size": 32},
         "experiment_name": "bs_32"},
        {"epochs": 100,
         "experiment_name": "epoch100"},
        {"lr": 0.001, "epochs": 50,
         "experiment_name": "lr_001"},
        {"lr": 0.001, "epochs": 50,
         "data": {"batch_size": 32},
         "experiment_name": "lr_001_bs_32"},
        {"lr": 0.01, "epochs": 50,
         "experiment_name": "lr_01"},
        {"lr": 0.01, "epochs": 50,
         "data": {"batch_size": 32},
         "experiment_name": "lr_01_bs_32"},
        {"lr": 0.00001, "epochs": 100,
         "experiment_name": "lr_00001"},
        {"model": {"layer_sizes": [29260, 300, 100, 10, 1]},
         "epochs": 50,
         "data": {"batch_size": 32},
         "experiment_name": "5layer_bs_32"},
        {"model": {"layer_sizes": [29260, 300, 100, 10, 1]},
         "epochs": 50,
         "data": {"batch_size": 16},
         "experiment_name": "5layer_bs_16"},
        {"model": {"layer_sizes": [29260, 10, 1]},
         "epochs": 50,
         "data": {"batch_size": 32},
         "experiment_name": "1_small_layer_bs_16"},
    ]
    return [apply_diff(BASE_CONFIG, diff) for diff in diffs]
