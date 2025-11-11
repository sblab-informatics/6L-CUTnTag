# ## import modules
import argparse
import random
import numpy as np
import pickle

import copy
import os

from scipy.stats import pearsonr, spearmanr

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, auc, average_precision_score, f1_score

from tqdm import tqdm

from sklearn.utils import shuffle

import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from dataset import HMCDataset

# customized modules
from model import enhancerHMC
# ================

def train(model, train_loader, lossfn, optimizer, schedular, device):
    model.train()
    total_loss = 0.
    preds, targets = [], []
    for batch, (enhancer, X, pad_mask, Y) in tqdm(enumerate(train_loader), desc='training', ncols=80):
    # for batch, (enhancer, X, pad_mask, Y) in enumerate(train_loader):
        X = X.to(device)
        pad_mask = pad_mask.to(device)
        Y = Y.to(device)

        pred, _ = model(X, pad_mask)

        loss = lossfn(pred, Y)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

        preds.extend(pred.detach().cpu().numpy().tolist())
        targets.extend(Y.cpu().numpy().tolist())
    
    schedular.step()

    preds_idxs = np.argmax(preds, axis=1)
    
    acc = accuracy_score(targets, preds_idxs)
    f1 = f1_score(targets, preds_idxs, average='macro')

    return total_loss / len(train_loader), acc, f1

# ================
def test(model, test_loader, lossfn, device):
    model.eval()
    preds, targets = [], []
    test_loss = 0.
    with torch.no_grad():
        for batch, (enhancer, X, pad_mask, Y) in tqdm(enumerate(test_loader), desc='testing', ncols=80):
        # for batch, (enhancer, X, pad_mask, Y) in enumerate(test_loader):
            X = X.to(device)
            pad_mask = pad_mask.to(device)
            Y = Y.to(device)

            pred, _ = model(X, pad_mask)

            loss = lossfn(pred, Y)
            test_loss += loss.item()

            preds.extend(pred.cpu().numpy().tolist())
            targets.extend(Y.cpu().numpy().tolist())
    
    preds_idxs = np.argmax(preds, axis=1)
    acc = accuracy_score(targets, preds_idxs)
    f1 = f1_score(targets, preds_idxs, average='macro')

    return test_loss / len(test_loader), acc, f1

# ===============

def load_data(args):

    fin = open('/home/xuan/Projects/MLfor6lCnT/data/WG_control_in_enhancers_modC.pkl', 'rb')
    # fin = open('/home/xuan/Projects/MLfor6lCnT/data/H3K4Me1_6L_in_enhancers_modC.pkl', 'rb')
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

    cpt_dir = f"checkpoints/checkpoints_with_random_seed{args.random_seed}"
    # pred_dir = f"predictions/predictions_with_random_seed{args.random_seed}"

    if not os.path.exists(cpt_dir):
        os.makedirs(cpt_dir)
    # if not os.path.exists(pred_dir):
    #     os.makedirs(pred_dir)

    tridxs, vaidxs, ttidxs = [], [], []
    for i in range(len(data[0])):
        en = data[0][i]
        chrom = en.split(':')[0]

        if chrom == 'chr1':
            ttidxs.append(i)

        elif chrom == 'chr10':
            vaidxs.append(i)
        
        else:
            tridxs.append(i)
    
    trdata = [data[j][tridxs] for j in range(len(data))]
    vadata = [data[j][vaidxs] for j in range(len(data))]
    ttdata = [data[j][ttidxs] for j in range(len(data))]

    trainset = HMCDataset(trdata)
    trainloader = DataLoader(trainset, batch_size=args.batch_size, num_workers=args.num_workers, shuffle=True, pin_memory=True)

    validset = HMCDataset(vadata)
    validloader = DataLoader(validset, batch_size=args.batch_size, num_workers=args.num_workers, shuffle=False, pin_memory=True)

    testset = HMCDataset(ttdata)
    testloader = DataLoader(testset, batch_size=args.batch_size, num_workers=args.num_workers, shuffle=False, pin_memory=True)

    enhancer_hmc = enhancerHMC(args).to(args.device)

    w1 = len(trdata[-1])/ np.sum(trdata[-1] == 0)
    w2 = len(trdata[-1])/ np.sum(trdata[-1] == 1)
    w3 = len(trdata[-1])/ np.sum(trdata[-1] == 2)
    Wcls = torch.tensor([w1, w2, w3], device=args.device).float()
    
    celoss = nn.CrossEntropyLoss(weight=Wcls)

    optimizer = torch.optim.Adam(enhancer_hmc.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    schedular = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, args.epochs, 0.7 * args.lr)

    # min_loss = 1e10
    max_acc = 0
    best_epoch = 0
    best_enhancer_hmc = copy.deepcopy(enhancer_hmc.state_dict())
    for epoch in range(args.epochs):
        trloss, tracc, trf1 = train(enhancer_hmc, trainloader, celoss, optimizer, schedular, args.device)
        valoss, vaacc, vaf1 = test(enhancer_hmc, validloader, celoss, args.device)
        print(f'Epoch:{epoch}; loss: {trloss:.4f}, {valoss:.4f}; acc: {tracc:.4f}, {vaacc:.4f}; f1: {trf1:.4f}, {vaf1:.4f}')
        print('='*80)
        
        # if ttloss < min_loss:
        if vaf1 > max_acc:
            # min_loss = ttloss
            max_acc = vaf1
            best_epoch = epoch
            best_enhancer_hmc = copy.deepcopy(enhancer_hmc.state_dict())
        else:
            if epoch - best_epoch >= args.patience:
                # print(f'=== Break at epoch {best_epoch + 1} ===')
                break
    
    enhancer_hmc.load_state_dict(best_enhancer_hmc)
        
    # cpt_path = f"{cpt_dir}/CnT6L_model_at_epoch{best_epoch+1}.cpt"
    cpt_path = f"{cpt_dir}/WG6L_model_at_epoch{best_epoch+1}.cpt"
    torch.save({
        'epoch': best_epoch + 1,
        'model_state_dict': enhancer_hmc.state_dict(),
        'acc': max_acc,
    }, cpt_path)
        
    # pred_path = f"{pred_dir}/model_at_epoch{best_epoch+1}_pred.pkl"
    ttloss, ttacc, ttf1 = test(enhancer_hmc, testloader, celoss, args.device)
    print(f'Epoch:{best_epoch+1}; loss: {ttloss:.4f}; acc: {ttacc:.4f}; f1: {ttf1:.4f}')

if __name__ == '__main__':

    main()

