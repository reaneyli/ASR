import SimpleITK as sitk
# from read_data import resize_data
import numpy as np
import os

def get_patch_D(imagepath,savepath):
    if not os.path.isdir(savepath):
        os.mkdir(savepath)
        print('create file directory finished！')
    filenames = os.listdir(imagepath)
    filenames.sort()
    COUNT = 0
    for filename in filenames:
        name, extension = os.path.splitext(filename)
        name1, extension1 = os.path.splitext(name)
        savepath1 = savepath+name+'/'
        if not os.path.isdir(savepath1):
            os.mkdir(savepath1)
            print('create file directory finished！')

        # print(name1)
        image = sitk.ReadImage(imagepath + filename)
        array = sitk.GetArrayFromImage(image)
        array = np.swapaxes(array, 0, 2)
        D, H, W = array.shape
        print(array.shape)  # (1533,145,100)
        ###################生成大小为（50，H,W）大小的patch
        # num = int(D / 80)
        D1 = 80
        width = 10
        for i in range(50):
            if i == 0:
                cropped = image[i * D1:(i + 1) * D1, :, :]
                if i>=1 and i<=9:
                    sitk.WriteImage(cropped, savepath1 + name1 + '_D0' + str(i) + '.mha')
                else:
                    sitk.WriteImage(cropped, savepath1 + name1  + '_D' + str(i) + '.mha')
                array_c = sitk.GetArrayFromImage(cropped)
                print('i:', i, array_c.shape)
            elif (i+1)*D1-i*width >= D:
                cropped = image[- D1:, :, :]
                if i>=1 and i<=9:
                    sitk.WriteImage(cropped, savepath1 + name1  + '_D0' + str(i) + '.mha')
                else:
                    sitk.WriteImage(cropped, savepath1 + name1  + '_D' + str(i) + '.mha')
                array_c = sitk.GetArrayFromImage(cropped)
                print('i:', i, array_c.shape)
                break
            else :
                cropped = image[i * D1 - i*width:(i + 1) * D1 - i*width, :, :]
                if i>=1 and i<=9:
                    sitk.WriteImage(cropped, savepath1 + name1 + '_D0' + str(i) + '.mha')
                else:
                    sitk.WriteImage(cropped, savepath1 + name1 + '_D' + str(i) + '.mha')
                array_c = sitk.GetArrayFromImage(cropped)
                print('i:', i, array_c.shape)
        COUNT += 1
        print('the', COUNT, 'images')

imagepath = 'G:/pytorch_celegons/data/test/testimgMHA/'
savepath = 'G:/pytorch_celegons/data/test/testimgpatch/'
get_patch_D(imagepath,savepath)
