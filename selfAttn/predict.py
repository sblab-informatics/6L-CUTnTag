# ## import modules
import argparse
import random
import numpy as np
import pickle

import copy
import os

from sklearn.metrics import accuracy_score, auc, average_precision_score, f1_score

from tqdm import tqdm

from sklearn.utils import shuffle

import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from dataset import HMCDataset

# customized modules
from model import TransHMC
# ================

def predict(model, pred_loader, opath, device):
    model.eval()
    enhancers_Wattns = {}
    with torch.no_grad():
        for batch, (enhancer, X, pad_mask, Y) in tqdm(enumerate(pred_loader), desc='testing', ncols=80):
        # for batch, (enhancer, X, pad_mask, Y) in enumerate(pred_loader):
            X = X.to(device)
            pad_mask = pad_mask.to(device)
            Y = Y.to(device)

            pred, Wattn = model(X, pad_mask)
            pad_mask = pad_mask.cpu().numpy()
            Wattn = Wattn.cpu().numpy()
            for i in range(len(enhancer)):
                hmc_non_empty = (pad_mask[i] == 0.).sum()
                w = Wattn[i][:hmc_non_empty, :hmc_non_empty]
                enhancers_Wattns[enhancer[i]] = w
    
    with open(opath, 'wb') as fout:
        pickle.dump(enhancers_Wattns, fout)

# ===============

def load_data(args):
    fin = open(args.input, 'rb')
    datadic = pickle.load(fin)
    fin.close()

    enhancers, counts, masks, labels = [], [], [], []
        
    for enhancer in tqdm(datadic.keys(), desc='loading data', ncols=80):
    # for enhancer in datadic.keys():
        
        modCs = datadic[enhancer] 

        if len(modCs) == 0:# or len(modCs) > args.max_hmc_num:
            continue

        mask = np.zeros(args.max_hmc_num)
        mask[len(modCs):] = 1.0
        masks.append(mask)

        enhancers.append(enhancer)
        
        l = enhancer.split(',')[1]
        if l == 'active':
            label = 0
        elif l == 'poised':
            label = 1
        else:
            label = 2

        labels.append(label)
        
        if len(modCs) > args.max_hmc_num:
            modCs = modCs[:args.max_hmc_num]
        
        if len(modCs) < args.max_hmc_num:
            l_pad = args.max_hmc_num - len(modCs)

            padCs = [[0.] * len(modCs[0])] * l_pad
            modCs.extend(padCs)

        counts.append(modCs)
    
    enhancers = np.array(enhancers)
    counts = np.array(counts)
    masks = np.stack(masks, axis=0)
    labels = np.array(labels)

    print("Shape of enhancers:", enhancers.shape)
    print("Shape of counts:", counts.shape)
    print("Shape of masks:", masks.shape)
    print("Shape of labels:", labels.shape)

    data = [enhancers, counts, masks, labels]

    data = shuffle(*data, random_state=args.random_seed)

    return data

def fix_random_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  
    np.random.seed(seed)  # Numpy module
    random.seed(seed)  # Python random module
    torch.backends.cudnn.deterministic = True
 
def get_args():
    parser = argparse.ArgumentParser(description='Attention-based model for enhancer classification from 6L-CnT')

    parser.add_argument('--input', type=str, default='', help='Path to the data dumped by pickle')
    parser.add_argument('--cpt', type=str, default='checkpoints/WG6L_model_at_epoch20.cpt', help='Path to the trained model')

    parser.add_argument('--device', type=int, default=0, help='GPU index')
    parser.add_argument('--num-workers', type=int, default=4)
    parser.add_argument('--batch-size', type=int, default=64)

    parser.add_argument('--hmc-d-model', type=int, default=32)
    parser.add_argument('--hmc-ff-dim', type=int, default=256)
    parser.add_argument('--hmc-act', type=str, default='elu')
    parser.add_argument('--hmc-n-head', type=int, default=4, help='Attention head number')
    parser.add_argument('--hmc-dropp', type=float, default=0.2)

    parser.add_argument('--hid-dim', type=int, default=256)
    parser.add_argument('--mlp-act', type=str, default='elu')

    parser.add_argument('--lr', type=float, default=1e-4, help='Learning rate')
    parser.add_argument('--weight-decay', type=float, default=0)
    parser.add_argument('--epochs', type=int, default=300)
    parser.add_argument('--patience', type=int, default=10)

    parser.add_argument('--max-hmc-num', type=int, default=100, help='maximum hmc sites considered')

    parser.add_argument('--random-seed', type=int, default=20250520, help='initial conditions for generating random variables')

    args = parser.parse_args()

    device = (f"cuda:{args.device}" if torch.cuda.is_available() else "cpu")
    args.device = device

    return args

# ================
def main():
    args = get_args()

    print(args)

    fix_random_seed(args.random_seed)
    data = load_data(args)

    pred_dir = f"predictions"

    if not os.path.exists(cpt_dir):
        os.makedirs(cpt_dir)
    if not os.path.exists(pred_dir):
        os.makedirs(pred_dir)

    allset = HMCDataset(data)
    allloader = DataLoader(allset, batch_size=args.batch_size, num_workers=args.num_workers, shuffle=False, pin_memory=True)
    state_dict = torch.load(f"{args.cpt}")['model_state_dict']

    transhmc = TransHMC(args).to(args.device)
    transhmc.load_state_dict(state_dict)

    opath = f"{pred_dir}/6LWG_predictions_attn.pkl"

    predict(transhmc, allloader, opath, args.device)

if __name__ == '__main__':

    main()

