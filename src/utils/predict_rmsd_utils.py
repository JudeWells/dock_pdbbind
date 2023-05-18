import os
import torch
import json

def basic_save_model(model, config, metrics={}):
    save_dir = os.path.join(config["checkpoint_dir"], config["experiment_name"])
    os.makedirs(save_dir, exist_ok=True)
    torch.save(model.state_dict(), os.path.join(*[save_dir, f'weights.pt']))
    with open(os.path.join(save_dir, 'config.json'), 'w') as f:
        json.dump(config, f)
    with open(os.path.join(save_dir, 'metrics.json'), 'w') as f:
        json.dump(metrics, f)


def save_val_best(metrics, learner, validation_metric='r2_score'):
    # TODO finish implementing this function
    current_metric = metrics.get(validation_metric, -float('inf'))
    if current_metric > learner.best_val_metrics.get(validation_metric, -float('inf')):
        best_val_save_dir = os.path.join(learner.checkpoint_dir, f'best_val')
        os.makedirs(best_val_save_dir, exist_ok=True)
        torch.save(learner.model.state_dict(), os.path.join(*[best_val_save_dir, f'weights.pt']))
        metrics_copy = {k:v for k,v in metrics.items()}
        metrics_copy['epoch'] = learner._epoch -1
        metrics_copy['validation_metric'] = validation_metric
        with open(os.path.join(best_val_save_dir, 'best_val_metrics.json'), 'w') as f:
            json.dump(metrics_copy, f)