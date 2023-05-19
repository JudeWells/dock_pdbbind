import sys

import torch
import wandb
from src.constructors import build_data_loaders
from src.models.docking_quality.ff import FeedForward
from src.configuration.base_config import BASE_CONFIG, first_run
from src.utils.predict_rmsd_utils import basic_save_model
from sklearn.metrics import r2_score

def evaluate_model(model, data_loader, data_loader_name):
    """
    Gets predictions for all data points in the dataloader and calculates
    The following metrics:
    - MSE
    - R2
    - MAE
    :param model:
    :param data_loader:
    :param data_loader_name:
    :return:
    """
    y_hat = []
    y = []
    with torch.no_grad():
        for i, (x, y_) in enumerate(data_loader):
            y_hat.append(model(x))
            y.append(y_)
        y_hat = torch.cat(y_hat)
        y = torch.cat(y)
        mse = float(torch.nn.functional.mse_loss(y_hat, y).detach())
        r2 = r2_score(y.detach().numpy(), y_hat.detach().numpy())
        mae = float(torch.nn.functional.l1_loss(y_hat, y).detach())
        metrics = {f"{data_loader_name}_mse": mse, f"{data_loader_name}_r2": r2, f"{data_loader_name}_mae": mae}
        print(f"{data_loader_name} MSE: {mse}, R2: {r2}, MAE: {mae}")
        if r2 > model.best_val:
            basic_save_model(model, config, metrics)
            model.best_val = r2
        wandb.log(metrics)
        return mse, r2, mae

if __name__ == "__main__":
    all_configs = first_run()
    try:
        task_index = int(sys.argv[1]) -1
        all_configs = [all_configs[task_index]]
    except:
        pass
    for i, config in enumerate(all_configs):
        train_loader, val_loader, test_loader = build_data_loaders(config.get("data", {}))
        model = FeedForward(layer_sizes=config["model"]["layer_sizes"], dropout_p=config["model"]["dropout_p"])

        print(model)
        optimizer = torch.optim.Adam(
            model.parameters(),
            lr=config.get("lr", 0.0001),
        )
        epochs = config.get("epochs", 15)
        run = wandb.init(
            config=config,
            entity="cath",
            project="rmse_dock",
            group=config.get("experiment_group", None),
            reinit=True,
            name=config["experiment_name"]
        )
        wandb.watch(model)
        for e in range(epochs):
            for i, (x, y) in enumerate(train_loader):
                optimizer.zero_grad()
                y_hat = model(x)
                loss = torch.nn.functional.mse_loss(y_hat, y)
                wandb.log({"batch_loss": loss.item()})
                loss.backward()
                optimizer.step()
            print(f"Epoch {e}, batch {i}, loss {loss.item()}")
            wandb.log({"epoch_loss": loss.item()})
            evaluate_model(model, val_loader, "val")
        evaluate_model(model, test_loader, "test")

