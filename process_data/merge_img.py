import numpy as np
from read_data import readimage,resize_data
import SimpleITK as sitk
import os
import math

#创建一个三维的全空数组，大小为（80，80，1600）
# 创建一张新图
#利用simpleITK将空数组存为图片，将这张图片分段赋值，最后重新保存

# （80，139，149）！！！

def get_imgshape_list(testimg_path):
    img_shape_list = []
    img_name_list = []
    filenames = os.listdir(testimg_path)
    filenames.sort()
    print(filenames)
    for i in range(len(filenames)):
        name, extension = os.path.splitext(filenames[i])
        name1, extension1 = os.path.splitext(name)
        image_array = readimage(testimg_path + filenames[i])
        print(image_array.shape)
        img_shape_list.append(image_array.shape)
        img_name_list.append(name1)

    return img_shape_list,img_name_list

def merge_D(filepath,savepath):
    filenames = os.listdir(filepath)
    filenames.sort()
    print(filenames)
    # new_array = np.zeros((1577, 149, 139))
    img_shape_list,img_name_list = get_imgshape_list(test_img_path)


    for i in range(len(filenames)):
        name, extension = os.path.splitext(filenames[i])
        new_array = np.zeros(img_shape_list[i])
        D, H, W = new_array.shape
        print(D, H, W)
        filenames1 = os.listdir(filepath+filenames[i])
        filenames1.sort()
        print(filenames1)
        for j in range(len(filenames1)):
            image = sitk.ReadImage(filepath + filenames[i]+'/'+ filenames1[j])
            array = sitk.GetArrayFromImage(image)

           # array = np.swapaxes(array, 0, 2)
         

            print(array.shape)
            

            width = 10
            D1 = 80
            if j == 0:
                new_array[:D1, :, :] = array
            elif (j + 1) * D1 - j * width > D:
                new_array[-D1:, :, :] = array
            else:
                new_array[j * D1 - j * width:(j + 1) * D1 - j * width] = array
        print(new_array.shape)
        new_array = np.swapaxes(new_array, 0, 2)
        new_image = sitk.GetImageFromArray(new_array)
        sitk.WriteImage(new_image, savepath + name + '.mha')

imgpath = '/home/qulab/wm/code4_unet/pred/'
test_img_path = '/home/qulab/wm/data/test1/testlabelMHA/'
savepath = '/home/qulab/wm/code4_unet/result/'
merge_D(imgpath,savepath)
