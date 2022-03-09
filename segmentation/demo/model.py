import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions.normal import Normal
import math


def make_one_hot(vol, mask):
    lens = len(mask)
    shape = np.array(vol.shape)
    shape[1] = lens
    shape = tuple(shape)
    result = torch.zeros(shape)

    for idx, label in enumerate(mask):
        tmp = vol == label
        result[:,idx] = tmp

    return result

def class2_diceloss(y_true, y_pred, smooth=1e-5):
    N = y_true.shape[1]
    y_pred_flat = y_pred.view(N,-1)
    y_true_flat = y_true.view(N,-1)
    intersection = y_pred_flat * y_true_flat
    loss = 2 * (intersection.sum() + smooth) / (y_pred_flat.sum() + y_true_flat.sum() + smooth)
    loss = 1 - loss.sum()/N
    return loss 


def dice_coefficient(y_true, y_pred, smooth=1e-5):
    y_true_d = np.sum(y_true*y_true)
    y_pred_d = np.sum(y_pred*y_pred)
    intersection = np.sum(y_true * y_pred)
    return (2. * intersection + smooth) / (y_true_d + y_pred_d + smooth)


def muldiceloss_1(y_true,y_pre):
    shape = y_true.shape
    total_loss = 0
    for i in range(shape[1]):
        dice_loss = class2_diceloss(y_true[:,i,:,:,:],y_pre[:,i,:,:,:])
        total_loss += dice_loss

    return total_loss

def dist_loss(y_true, y_pred, smooth=1e-5):
    loss = (y_true - y_pred) * (y_true - y_pred)/(y_true.shape[0] *y_true.shape[1] *y_true.shape[2] * y_true.shape[3]  * y_true.shape[4])
    return loss.sum()

def dist_loss_mul(y_true, y_pred):
    shape = y_true.shape

    total_loss = 0
    for i in range(shape[1]):
        loss = (y_true[:,i,:,:,:] - y_pred[:,i,:,:,:]) * (y_true[:,i,:,:,:] - y_pred[:,i,:,:,:])
        total_loss += loss.sum()
    total_loss = total_loss/(y_true.shape[0] *y_true.shape[1] *y_true.shape[2] * y_true.shape[3]  * y_true.shape[4])
   
    return total_loss

def color_loss(y_true,y_pred,y):

    similarity = torch.cosine_similarity(y_true, y_pred, dim=1)
    similarity[y==0]=1
    similarity=1-similarity
    return torch.sum(similarity)/(y_true.shape[0]  *y_true.shape[2] * y_true.shape[3]  * y_true.shape[4])
