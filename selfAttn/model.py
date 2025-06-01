import torch
import torch.nn as nn

import torchsnooper

class AttnBlock(nn.Module):
    def __init__(self, d_model, n_head, ff_dim, act_fn, dropp):
        super().__init__()

        self.attn = nn.MultiheadAttention(d_model, n_head, batch_first=True)
        
        self.ffn = nn.Sequential(
            nn.Linear(d_model, ff_dim),
            act_fn(),
            nn.Linear(ff_dim, d_model)
        )

        self.ln1 = nn.LayerNorm(d_model)
        self.ln2 = nn.LayerNorm(d_model)

        self.drop1 = nn.Dropout(dropp)
        self.drop2 = nn.Dropout(dropp)
    
    def forward(self, x, pad_mask=None):
        x_att, w_att = self.attn(x, x, x, key_padding_mask=pad_mask)
        x = self.ln1(x + self.drop1(x_att))

        x_ffn = self.ffn(x)
        x = self.ln2(x + self.drop2(x_ffn))

        return x, w_att

class TransHMC(nn.Module):
    def __init__(self, args):
        super().__init__()

        act_dic = {
            'elu': nn.ELU,
            'silu': nn.SiLU,
            'relu': nn.ReLU,
            'gelu': nn.GELU,
            'prelu': nn.PReLU,
            'lrelu': nn.LeakyReLU
        }

        hmc_attn_act = act_dic[args.hmc_act]
        mlp_act = act_dic[args.mlp_act]

        self.hmc_proj = nn.Linear(4, args.hmc_d_model)
        
        self.hmc_block = AttnBlock(args.hmc_d_model, args.hmc_n_head, args.hmc_ff_dim, hmc_attn_act, args.hmc_dropp)

        self.flat = nn.Flatten()

        self.mlp = nn.Sequential(
            nn.Linear(args.hmc_d_model * args.max_hmc_num, args.hid_dim),
            mlp_act(),
            nn.Linear(args.hid_dim, args.hid_dim),
            mlp_act(),
            nn.Linear(args.hid_dim, 3),
            # nn.Softmax(dim=-1)
        )

    # @torchsnooper.snoop()
    def forward(self, x_hmc, hmc_pad_mask):
        x = self.hmc_proj(x_hmc)

        x, attn_w = self.hmc_block(x, hmc_pad_mask)

        x = self.flat(x)
        
        x = self.mlp(x)

        return x, attn_w
