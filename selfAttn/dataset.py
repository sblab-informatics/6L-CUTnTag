import torch
from torch.utils.data import Dataset
import numpy as np

class HMCDataset(Dataset):
    """
    members:
        data: [enhancer, hmc_counts, hmc_masks, cls]
    """

    def __init__(self, data):
        super().__init__()

        enhancer, X, pad_masks, Y = data

        self.enhancer = enhancer
        self.X = torch.from_numpy(X).float()
        self.pad_masks = torch.from_numpy(pad_masks).float()
        self.Y = torch.from_numpy(Y).long()

    def __len__(self):
        return len(self.Y)
    
    def __getitem__(self, idx):
        return self.enhancer[idx], self.X[idx], self.pad_masks[idx], self.Y[idx]