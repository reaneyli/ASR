import os

from scipy import ndimage
import SimpleITK as sitk
from read_data import readimage


datapath = 'G:/pytorch_celegons/data/train/train_label_binary/'
savepath = 'G:/pytorch_celegons/data/train/train_label_edt/'

def save_edt_img(datapath,savepath):
    filenames = os.listdir(datapath)
    filenames.sort()
    print(filenames)
    for i in range(len(filenames)):
        img = datapath + filenames[i]
        input = readimage(img)
        distance_array = edt_trans(input)
        path, fullname = os.path.split(img)
        name, extension = os.path.splitext(fullname)
        name1, extension1 = os.path.splitext(name)

        # print('name, extension ',name,extension)
        # print('name1, extension1',name1, extension1)
        # print('name2, extension2', name2, extension2)
        image = sitk.GetImageFromArray(distance_array)

        sitk.WriteImage(image, savepath + name1 + '_edt' + '.mha')
        print('finish!', i)


def BianryImage(ImageArray):
    ImageArray[ImageArray == 0] = 0
    ImageArray[ImageArray > 0] = 255

    return ImageArray

def edt_trans(input):
    distance_array = ndimage.morphology.distance_transform_edt(input, sampling=None, return_distances=True,
                                                               return_indices=False, distances=None, indices=None)
    return distance_array

def get_center_mask(distance_array):
    # imgarray = readimage(img)
    # distance_array[distance_array==0]=0
    distance_array[distance_array >1]=255
    distance_array[distance_array<1] = 0
    return distance_array

def get_inside_mask(distance_array):
    distance_array[distance_array == 0] = 0
    distance_array[distance_array > 1] = 0
    distance_array[distance_array != 0] = 255

    return distance_array
# img = 'G:/pytorch_celegons/data/train/labelMHA/aha-1_XIL029_20131021_0004_mask.raw.mha'
# input = readimage(img)
# distance_array = edt_trans(input)
# distance_array1 = edt_trans(input)
# center_mask = get_center_mask(distance_array)
# inside_mask = get_inside_mask(distance_array1)
# image = sitk.GetImageFromArray(center_mask)
# image1 = sitk.GetImageFromArray(inside_mask)
# sitk.WriteImage(image, savepath  + 'center' + '.mha')
# sitk.WriteImage(image1, savepath  + 'inside' + '.mha')

