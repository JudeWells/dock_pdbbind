from torch import nn


class FeedForward(nn.Module):
    def __init__(self, layer_sizes, dropout_p=0.15):
        super().__init__()
        layers = []
        for i in range(len(layer_sizes) - 1):
            input_size = layer_sizes[i]
            output_size = layer_sizes[i + 1]
            layers.append(nn.Linear(input_size, output_size))
            if i < len(layer_sizes) - 2:  # no batchnorm/activation/dropout on last layer
                layers.append(nn.BatchNorm1d(output_size))
                layers.append(nn.ELU())
                layers.append(nn.Dropout(p=dropout_p))
        self.layers = nn.Sequential(*layers)

    def forward(self, x):
        return self.layers(x)
