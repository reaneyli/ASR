import SimpleITK as sitk
import torchvision.transforms as tfs
import numpy as np
import torch
import os
from PIL import Image
from scipy.ndimage import zoom
from torch.utils import data

class C_Dataset(data.Dataset):
    def __init__(self,datapath):

        imgs =[]
        filenames = os.listdir(datapath)
        filenames.sort()
        print(filenames)
        filenames_img = os.listdir(datapath+filenames[0])
        filenames_img.sort()
        filenames_label = os.listdir(datapath+filenames[1])
        filenames_label.sort()
        for i in range(len(filenames_img)):
            img = datapath + filenames[0] + '/' + filenames_img[i]
            label = datapath +filenames[1] + '/' + filenames_label[i]
            imgs.append([img,label])



        self.imgs = imgs

    def __getitem__(self, index):
        x_path,y_path = self.imgs[index]
        img_x = readimage(x_path)
        img_x = resize_data(img_x,(128, 128, 80))
        img_x = np.reshape(img_x, (1,) + img_x.shape)



        img_y = sitk.ReadImage(y_path)
        img_y = sitk.GetArrayFromImage(img_y)
        img_y = resize_data_label(img_y,(3,128, 128, 80))


        x_transform = torch.Tensor(img_x)
        y_transform = torch.Tensor(img_y)
        return x_transform,y_transform

    def __len__(self):
        return len(self.imgs)



def readimage(imagepath):
    image = sitk.ReadImage(imagepath)
    array = sitk.GetArrayFromImage(image)
    array = np.swapaxes(array,0,2)
    return array

def resize_data(data,input_shape):
    (input_D, input_H, input_W) = input_shape
    [depth, height, width] = data.shape
    scale = [input_D * 1.0 / depth, input_H * 1.0 / height, input_W * 1.0 / width]
    data = zoom(data, scale, order=3)
    return data

def resize_data_label(data,input_shape):
    (input_C,input_D, input_H, input_W) = input_shape
    [c,depth, height, width] = data.shape
    scale = [input_C* 1.0 / c,input_D * 1.0 / depth, input_H * 1.0 / height, input_W * 1.0 / width]
    data = zoom(data, scale, order=0)
    return data


def BianryImage(ImageArray):
    ImageArray[ImageArray == 0] = 0
    ImageArray[ImageArray > 0] = 1
    return ImageArray

