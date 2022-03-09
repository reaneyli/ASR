#!/usr/bin/env python
# -*- coding: utf-8 -*-
import SimpleITK as sitk
from skimage.morphology import watershed
import numpy as np
from skimage.measure import label
from skimage.morphology import reconstruction, dilation, erosion, disk, diamond, square
from skimage import img_as_ubyte
import time
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.metrics import jaccard_score, f1_score
from sklearn.metrics import recall_score, precision_score, confusion_matrix
from skimage.morphology import erosion, disk
from os.path import join
from skimage.io import imsave, imread
import os
from compute_score import scores,scores_AJI
from read_data import readimage
from fill_hole import fill



def dice_coefficient(y_true, y_pred, smooth=1e-5):
    y_true_d = np.sum(y_true*y_true)
    y_pred_d = np.sum(y_pred*y_pred)
    intersection = np.sum(y_true * y_pred)
    return (2. * intersection + smooth) / (y_true_d + y_pred_d + smooth)

def PrepareProb(img, convertuint8=True, inverse=True):
    """
    Prepares the prob image for post-processing, it can convert from
    float -> to uint8 and it can inverse it if needed.
    """
    if convertuint8:
        img = img_as_ubyte(img)
    if inverse:
        img = 255 - img
    return img


def HreconstructionErosion(prob_img, h):
    """
    Performs a H minimma reconstruction via an erosion method.
    """

    def making_top_mask(x, lamb=h):
       return min(255, x + lamb)

    f = np.vectorize(making_top_mask)
    shift_prob_img = f(prob_img)  #返回prob_img + h和255之间的最小值
    # shift_prob_img = sitk.GetImageFromArray(shift_prob_img)
    # sitk.WriteImage(shift_prob_img, savepath + 'shift_prob_img' + '.mha')
    seed = shift_prob_img
    mask = prob_img
    recons = reconstruction(
      seed, mask, method='erosion').astype(np.dtype('ubyte'))   #返回图像的形态学重构（在这里通过腐蚀进行重建）
    return recons


def find_maxima(img, convertuint8=False, inverse=False, mask=None):  #找到图像的所有局部最大值
    """
    Finds all local maxima from 2D image.
    """
    img = PrepareProb(img, convertuint8=convertuint8, inverse=inverse)
    recons = HreconstructionErosion(img, 1)
    if mask is None:
        return recons - img
    else:
        res = recons - img
        res[mask==0] = 0
        return res

def GetContours(img):
    """
    Returns only the contours of the image.
    The image has to be a binary image 
    """
    img[img > 0] = 1
    return dilation(img, disk(2)) - erosion(img, disk(2))


def generate_wsl(ws):
    """
    Generates watershed line that correspond to areas of
    touching objects.
    """
    print('ws.shape',ws.shape)
    se = square(3)
    se1 = [[[1, 1, 1], [1, 1, 1], [1, 1, 1]], [[1, 1, 1], [1, 1, 1], [1, 1, 1]], [[1, 1, 1], [1, 1, 1], [1, 1, 1]]]
    se1 = np.array(se1)
    ero = ws.copy()
    ero[ero == 0] = ero.max() + 1
    ero = erosion(ero, se1)
    ero[ws == 0] = 0

    grad = dilation(ws, se1) - ero
    grad[ws == 0] = 0
    grad[grad > 0] = 255
    grad = grad.astype(np.uint8)
    return grad

def DynamicWatershedAlias(p_img, lamb, p_thresh = 0.5):
    """
    Applies our dynamic watershed to 2D prob/dist image.
    """
    b_img = (p_img > p_thresh) + 0
    Probs_inv = PrepareProb(p_img)

    Hrecons = HreconstructionErosion(Probs_inv, lamb)

    markers_Probs_inv = find_maxima(Hrecons, mask = b_img)

    markers_Probs_inv = label(markers_Probs_inv)
    markers_Probs_inv1 = markers_Probs_inv.astype(np.float)

    ws_labels = watershed(Hrecons, markers_Probs_inv, mask=b_img)
    arrange_label = ArrangeLabel(ws_labels)

    arrange_label1 = arrange_label.astype(np.float)

    wsl = generate_wsl(arrange_label)

    arrange_label[wsl > 0] = 0
    arrange_label2 = arrange_label.astype(np.float)

    return arrange_label

def ArrangeLabel(mat):
    """
    Arrange label image as to effectively put background to 0.
    """
    val, counts = np.unique(mat, return_counts=True)
    background_val = val[np.argmax(counts)]
    mat = label(mat, background = background_val)
    if np.min(mat) < 0:
        mat += np.min(mat)
        mat = ArrangeLabel(mat)
    return mat


def PostProcess(prob_image, param=7, thresh = 0.5):
    """
    Perform DynamicWatershedAlias with some default parameters.
    """
    segmentation_mask = DynamicWatershedAlias(prob_image, param, thresh)
    return segmentation_mask

def ComputeMetrics(prob, batch_labels, p1, p2, rgb=None, save_path=None, ind=0):  #这里p1，p2设置的分别为7和0.5
    """
    Computes all metrics between probability map and corresponding label.
    """
    GT = label(batch_labels.copy())
    PRED = prob
    PRED = label(PRED)
    #若处理的是需要后处理的图片，则把以下三行取消注释
    # PRED = PostProcess(prob, p1, p2)
    # PRED = fill(PRED)
    # PRED = label((prob > 0.5).astype('uint8'))
    lbl = GT.copy()
    pred = PRED.copy()
    lbl[lbl > 0] = 1
    pred[pred > 0] = 1
    dice = dice_coefficient(lbl,pred)
    l, p = lbl.flatten(), pred.flatten()
    acc = accuracy_score(l, p)
    roc = roc_auc_score(l, p)
    jac = jaccard_score(l, p)
    f1 = f1_score(l, p)
    recall = recall_score(l, p)
    precision = precision_score(l, p)


    return acc, roc, jac, recall, precision, f1,PRED,dice


def metrics(labelpath,prepath,p1,p2):
    label1 = sitk.ReadImage(labelpath)
    ground_truth = sitk.GetArrayFromImage(label1)
    print('ground_truth.shape:',ground_truth.shape)
    GT = label(ground_truth.copy())
    GT = GT.astype(np.float32)
    print('GT.max:',GT.max())

    pre = sitk.ReadImage(prepath)
    prediction = sitk.GetArrayFromImage(pre)
    #prediction = np.swapaxes(prediction,0,2)
    print('prediction.shape:', prediction.shape)

    # 如果是测试需要后处理的图片，则把以下三行取消注释
    # prediction[prediction<0.1]=0
    # prediction = prediction / prediction.max()
    # print('prediction.max:',prediction.max())

    acc, roc, jac, recall, precision, f1,PRED ,dice= ComputeMetrics(prediction, ground_truth, p1, p2)
    cell_num = PRED.max()
    print(acc, roc, jac, recall, precision, f1,cell_num)
    return PRED,acc, roc, jac, recall, precision,dice, f1,cell_num


prepath = 'F:/segmentation_mha/'
labelpath = 'G:/pytorch_celegons/data/test/testlabelMHA/'
savepath = 'G:/pytorch_celegons/try/5/'


filter = [[[1,1,1],[1,1,1],[1,1,1]],[[1,1,1],[1,1,1],[1,1,1]],[[1,1,1],[1,1,1],[1,1,1]]]
filter = np.array(filter)

# p1 = 22
p2  = 0.0


vis_array = []
for p1 in range(16,17):
    sum_acc = 0
    sum_roc = 0
    sum_recall = 0
    sum_precision = 0
    sum_dice = 0
    sum_f1 = 0
    sum_cell_num = 0
    count = 0

    sum_correct = 0
    sum_error = 0
    sum_score = 0
    sum_AJI = 0
    sum_pre_single = 0

    print("p1=",p1)
    single_array = []
    label_filenames = os.listdir(labelpath)
    label_filenames.sort()
    print(label_filenames)

    pre_filenames = os.listdir(prepath)
    pre_filenames.sort()
    print(pre_filenames)
    for i in range(len(pre_filenames)):
        count += 1
        PRED, acc, roc, jac, recall, precision, dice,f1, cell_num = metrics(labelpath + label_filenames[i],
                                                                 prepath + pre_filenames[i], p1, p2)

        label_array = readimage(labelpath+ label_filenames[i])
        print('label_array.shape',label_array.shape)
        # PRED = np.swapaxes(PRED,0,2)
        print('PRED.shape', PRED.shape)
        # **************************************
        if not os.path.isdir(savepath):
            os.mkdir(savepath)

        print('111111111111111111111111111:', PRED.shape, PRED.dtype)
        PRED = PRED.astype(np.uint16)
        print('111111111111111111111111111:', PRED.shape, PRED.dtype)
        # Img = sitk.GetImageFromArray(PRED, isVector=False)
        # sitk.WriteImage(Img, savepath + pre_filenames[i])


        correct, error, score,pre_single,sum_jiaoji_count,sum_bingji_count,cell_array = scores_AJI(label_array, PRED)
        AJI = sum_jiaoji_count/(sum_bingji_count+pre_single)
        sum_acc += acc
        sum_roc += roc
        sum_recall += recall
        sum_precision += precision
        sum_dice += dice
        sum_f1 += f1
        sum_cell_num += cell_num

        sum_correct += correct
        sum_error += error
        sum_score += score
        sum_AJI += AJI
        sum_pre_single += pre_single

        print(p1,acc, roc, jac, recall, precision, f1,cell_num,correct, error, score,AJI)
        print('...........................', i,'model')
    aver_acc = sum_acc / count
    aver_roc = sum_roc / count
    aver_recall = sum_recall / count
    aver_precision = sum_precision / count
    aver_dice = sum_dice/count
    aver_f1 = sum_f1 / count
    aver_cell_num = sum_cell_num / count

    aver_correct = sum_correct/count
    aver_error = sum_error/count
    aver_score = sum_score/count
    aver_AJI = sum_AJI/count
    aver_pre_single = sum_pre_single/count

    single_array.append([p1,aver_acc,aver_roc,aver_recall,aver_precision,aver_dice,aver_f1,aver_cell_num,aver_correct,aver_error,aver_score,aver_AJI,aver_pre_single])
    vis_array.append([p1, aver_acc, aver_roc, aver_recall, aver_precision,aver_dice, aver_f1, aver_cell_num,aver_correct,aver_error,aver_score,aver_AJI,aver_pre_single])
    # np.savetxt('merge_patch_8_1_epoch_much_' + str(p1)+ '.txt',vis_array,fmt='%.5f')
    # np.savetxt('cell_array_unet' + str(p1)+ '.txt',cell_array,fmt='%.5f')
    print('merge_result_7_1 p1:', p1)
    print('aver_acc:', aver_acc)
    print('aver_roc:', aver_roc)
    print('aver_recall:', aver_recall)
    print('aver_precision:', aver_precision)
    print('aver_dice:', aver_dice)
    print('aver_f1:', aver_f1)
    print('aver_cell_num:', aver_cell_num)
    

print(vis_array)
np.savetxt('vis_array.txt',vis_array,fmt='%.5f')


