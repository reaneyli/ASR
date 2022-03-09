import SimpleITK as sitk
import numpy as np
import torch
import os
import tables
import h5py

from scipy.ndimage import zoom


def readimage(imagepath):
    image = sitk.ReadImage(imagepath)
    array = sitk.GetArrayFromImage(image)
    array = np.swapaxes(array,0,2)
    return array


def read_h5_file(filepath, readwrite='r'):
    return tables.open_file(filepath, readwrite)

def resizedata(data,input_shape):
    # (input_D, input_H, input_W,input_C) = input_shape
    # [depth, height, width, channel] = data.shape
    # scale = [input_D * 1.0 / depth, input_H * 1.0 / height, input_W * 1.0 / width, input_C * 1.0 / channel]
    # data = zoom(data, scale, order=0)

    (input_D, input_H, input_W) = input_shape
    [depth, height, width] = data.shape
    # print(data.max())
    scale = [input_D * 1.0 / depth, input_H * 1.0 / height, input_W * 1.0 / width]
    data = zoom(data, scale, order=0)
    # print(data.max())
    return data



# datapath ='G:/celegans/data/train/'
# filenames = os.listdir(datapath)
# filenames.sort()
# filenames_img = os.listdir(datapath+filenames[1])
# filenames_img.sort()
# filenames_label = os.listdir(datapath+filenames[3])
# filenames_label.sort()
# for i in range(len(filenames_img)):
#     img = datapath + filenames[1] + '/' + filenames_img[i]
#     image =readimage(img)
#     print(im