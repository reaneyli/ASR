import os
from argparse import ArgumentParser
import numpy as np
import torch
import UNET
from read_data import readimage,resize_data
import SimpleITK as sitk

parser = ArgumentParser()
parser.add_argument("--modelpath", type=str,
                    dest="modelpath", default='/home/qulab/wm/code4_unet/MODEL/unet3d.pth',
                    help="frequency of saving models")
parser.add_argument("--savepath", type=str,
                    dest="savepath", default='./pred/',
                    help="path for saving images")
parser.add_argument("--start_channel", type=int,
                    dest="start_channel", default=16,
                    help="number of start channels")
args = parser.parse_args()

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if torch.cuda.is_available():
    print('cuda is available')
savepath = args.savepath

if not os.path.isdir(savepath):
    os.mkdir(savepath)

path = '/home/qulab/wm/data/test1/testimgpatch/'

def save_img(img,shape, savename):
    img = resize_data(img, shape)
    print(img.shape)
    Img = sitk.GetImageFromArray(img, isVector=False)
    sitk.WriteImage(Img, savename)


def visualize():
    model = UNET.unet3d(1, 3)
    model.to(device)
   
    model.load_state_dict(torch.load(args.modelpath),False)
    model.eval()

    imgs =[]
    filenames = os.listdir(path)
    filenames.sort()
    # print(filenames)

    for i in range(len(filenames)):
        filenames1 = os.listdir(path + filenames[i])
        filenames1.sort()
        # print(filenames1)
        savepath1 = savepath + filenames[i] + '/'
        if not os.path.isdir(savepath1):
            os.mkdir(savepath1)
            print("finished create!!!")
        for j in range(len(filenames1)):
            print('------------------------------',i,j)
            x = readimage(path + filenames[i]+'/' + filenames1[j])
            x = np.swapaxes(x, 0, 2)
            shape = x.shape
            print(x.shape)
            img_y = resize_data(x, (128, 128, 80))
            img_y = np.reshape(img_y, (1,) + (1,) + img_y.shape)
            img_y = torch.tensor(img_y, dtype=torch.float)
            print(img_y.shape)
            img_y = img_y.to(device)
            y_pre = model(img_y)
            print('y_pre shape :', y_pre.shape)
            print('y_pre.max= :', y_pre.max())
            y_pre = y_pre * y_pre
            print('1 :', y_pre.shape)
            y_pre = torch.sqrt(y_pre[:, 0, :, :, :] + y_pre[:, 1, :, :, :] + y_pre[:, 2, :, :, :])
            print('2 :', y_pre.shape)
            print('2  :', y_pre.max())
            print('2  :', y_pre.min())
            y_pre[y_pre < 0.4] = 0
            y_pre = torch.squeeze(y_pre)
            print('y_pre shape :', y_pre.shape)

            y_pre_save = y_pre.permute(2, 1, 0).detach().cpu().numpy()[:, :, :]
            print('y_pre_save shape:', y_pre_save.shape)
            y_pre_save = resize_data(y_pre_save, (shape[2], shape[1], shape[0]))
            print('y_pre_save shape:', y_pre_save.shape)
            Img = sitk.GetImageFromArray(y_pre_save, isVector=False)
            sitk.WriteImage(Img, savepath1 + filenames1[j])

visualize()

