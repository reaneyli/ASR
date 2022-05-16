import SimpleITK as sitk
# from read_data import resize_data
import numpy as np
import os
import struct
from scipy.ndimage import zoom

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

def load_v3d_raw_img_file(filename):
    im = {}
    try:
        f_obj = open(filename, 'rb')
    except FileNotFoundError:
        print("ERROR: Failed in reading [" + filename + "], Exit.")
        f_obj.close()
        return im
    else:
        # read image header - formatkey(24bytes)
        len_formatkey = len('raw_image_stack_by_hpeng')
        formatkey = f_obj.read(len_formatkey)
        formatkey = struct.unpack(str(len_formatkey) + 's', formatkey)
        if formatkey[0] != b'raw_image_stack_by_hpeng':
            print("ERROR: File unrecognized (not raw, v3draw) or corrupted.")
            f_obj.close()
            return im

        # read image header - endianCode(1byte)
        endiancode = f_obj.read(1)
        endiancode = struct.unpack('c', endiancode)  # 'c' = char
        endiancode = endiancode[0]
        if endiancode != b'B' and endiancode != b'L':
            print("ERROR: Only supports big- or little- endian,"
                  " but not other format. Check your data endian.")
            f_obj.close()
            return im

        # read image header - datatype(2bytes)
        datatype = f_obj.read(2)
        if endiancode == b'L':
            datatype = struct.unpack('<h', datatype)  # 'h' = short
        else:
            datatype = struct.unpack('>h', datatype)  # 'h' = short
        datatype = datatype[0]
        if datatype < 1 or datatype > 4:
            print("ERROR: Unrecognized data type code [%d]. "
                  "The file type is incorrect or this code is not supported." % (datatype))
            f_obj.close()
            return im

        # read image header - size(4*4bytes)
        size = f_obj.read(4 * 4)
        if endiancode == b'L':
            size = struct.unpack('<4l', size)  # 'l' = long
        else:
            size = struct.unpack('>4l', size)  # 'l' = long
        # print(size)

        # read image data
        npixels = size[0] * size[1] * size[2] * size[3]
        im_data = f_obj.read()
        if datatype == 1:
            im_data = np.frombuffer(im_data, np.uint8)
        elif datatype == 2:
            im_data = np.frombuffer(im_data, np.uint16)
        else:
            im_data = np.frombuffer(im_data, np.float32)
        if len(im_data) != npixels:
            print("ERROR: Read image data size != image size. Check your data.")
            f_obj.close()
            return im

        im_data = im_data.reshape((size[3], size[2], size[1], size[0]))
        # print(im_data.shape)
        im_data = np.moveaxis(im_data, 0, -1)
        # print(im_data.shape)
        im_data = np.moveaxis(im_data, 0, -2)
        # print(im_data.shape)
    f_obj.close()

    im['endian'] = endiancode
    im['datatype'] = datatype
    im['size'] = im_data.shape
    im['data'] = im_data
    return im
def get_patch_D(imagepath,savepath,savepath_r,data_form):
    if not os.path.isdir(savepath):
        os.mkdir(savepath)
        print('create file directory finished！')
    filenames = os.listdir(imagepath)
    filenames.sort()
    COUNT = 0
    for filename in filenames:
        name, extension = os.path.splitext(filename)
        name1, extension1 = os.path.splitext(name)
        savepath1 = savepath+name1+'/'
        savepath2 = savepath_r + name1 + '/'
        if not (os.path.isdir(savepath1) or os.path.isdir(savepath2)):
            os.mkdir(savepath1)
            os.mkdir(savepath2)
            print('create file directory finished！')

        # print(name1)
        if (extension==".mha"):
            image = sitk.ReadImage(imagepath + filename)
            array = sitk.GetArrayFromImage(image)
            #array = np.swapaxes(array, 0, 2)
        elif (extension==".v3draw"):
            image = load_v3d_raw_img_file(imagepath + filename)
            array = np.squeeze(image['data'])
            array = np.swapaxes(array, 0, 1)

        C,D, H, W = array.shape
        print(array.shape)  # (1533,145,100)
        ###################生成大小为（50，H,W）大小的patch
        # num = int(D / 80)
        D1 = 80
        width = 10
        for i in range(50):
            if i == 0:

                if (extension == ".mha"):
                    cropped = array[:,i * D1:(i + 1) * D1, :, :]
                    cropped = np.swapaxes(cropped, 1, 3)
                    img_y = resize_data_label(cropped, (3, 128, 128, 80))
                    img_y = sitk.GetImageFromArray(img_y, isVector=False)

                elif (extension == ".v3draw"):
                    if(data_form=="image"):
                        cropped = array[i * D1:(i + 1) * D1, :, :]
                    elif(data_form=="label"):
                        cropped1 = array[i * D1:(i + 1) * D1, :, :]
                        cropped2 =array[int(D / 3) + i * D1:int(D / 3) +(i + 1) * D1, :, :]
                        cropped3 = array[int(D / 3*2) + i * D1:int(D / 3*2) +(i + 1) * D1, :, :]
                        cropped = np.append(cropped1, cropped2, axis=0)
                        cropped = np.append(cropped, cropped3, axis=0)
                        cropped = np.swapaxes(cropped, 0, 2)
                cropped = sitk.GetImageFromArray(cropped, isVector=False)
                if i>=1 and i<=9:
                    sitk.WriteImage(cropped, savepath1 + name1 + '_D0' + str(i) + '.mha')
                    sitk.WriteImage(img_y, savepath2 + name1 + '_D0' + str(i) + '.mha')
                else:
                    sitk.WriteImage(cropped, savepath1 + name1  + '_D' + str(i) + '.mha')
                    sitk.WriteImage(img_y, savepath2 + name1 + '_D' + str(i) + '.mha')
                array_c = sitk.GetArrayFromImage(cropped)
                print('i:', i, array_c.shape)
            elif (i+1)*D1-i*width >= D:

                if (extension == ".mha"):
                    cropped = array[:,- D1:, :, :]
                    cropped = np.swapaxes(cropped, 1, 3)
                    img_y = resize_data_label(cropped, (3, 128, 128, 80))
                    img_y = sitk.GetImageFromArray(img_y, isVector=False)
                elif (extension == ".v3draw"):
                    if (data_form == "image"):
                        cropped = array[- D1:, :, :]
                    elif (data_form == "label"):
                        cropped1 = array[- D1:, :, :]
                        cropped2 = array[int(D / 3*2) - D1:int(D / 3*2), :, :]
                        cropped3 = array[int(D ) - D1:int(D ), :, :]
                        cropped = np.append(cropped1, cropped2, axis=0)
                        cropped = np.append(cropped, cropped3, axis=0)
                        cropped = np.swapaxes(cropped, 0, 2)
                cropped = sitk.GetImageFromArray(cropped, isVector=False)
                if i>=1 and i<=9:
                    sitk.WriteImage(cropped, savepath1 + name1  + '_D0' + str(i) + '.mha')
                    sitk.WriteImage(img_y, savepath2 + name1 + '_D0' + str(i) + '.mha')
                else:
                    sitk.WriteImage(cropped, savepath1 + name1  + '_D' + str(i) + '.mha')
                    sitk.WriteImage(img_y, savepath2 + name1 + '_D' + str(i) + '.mha')
                array_c = sitk.GetArrayFromImage(cropped)
                print('i:', i, array_c.shape)
                break
            else :

                if (extension == ".mha"):
                    cropped = array[:,i * D1 - i * width:(i + 1) * D1 - i * width, :, :]
                    cropped = np.swapaxes(cropped, 1, 3)
                    img_y = resize_data_label(cropped, (3, 128, 128, 80))
                    img_y = sitk.GetImageFromArray(img_y, isVector=False)
                elif (extension == ".v3draw"):
                    if(data_form=="image"):
                        cropped = array[i * D1 - i * width:(i + 1) * D1 - i * width, :, :]
                    elif(data_form=="label"):
                        cropped1 = array[i * D1 - i * width:(i + 1) * D1 - i * width, :, :]
                        cropped2 =array[int(D / 3) + i * D1 - i * width:int(D / 3) +(i + 1) * D1 - i * width, :, :]
                        cropped3 = array[int(D / 3*2) + i * D1 - i * width:int(D / 3*2) +(i + 1) * D1 - i * width, :, :]
                        cropped = np.append(cropped1, cropped2, axis=0)
                        cropped = np.append(cropped, cropped3, axis=0)
                        cropped = np.swapaxes(cropped, 0, 2)
                cropped = sitk.GetImageFromArray(cropped, isVector=False)

                if i>=1 and i<=9:
                    sitk.WriteImage(cropped, savepath1 + name1 + '_D0' + str(i) + '.mha')
                    sitk.WriteImage(img_y, savepath2 + name1 + '_D0' + str(i) + '.mha')
                else:
                    sitk.WriteImage(cropped, savepath1 + name1 + '_D' + str(i) + '.mha')
                    sitk.WriteImage(img_y, savepath2 + name1 + '_D' + str(i) + '.mha')
                array_c = sitk.GetArrayFromImage(cropped)
                print('i:', i, array_c.shape)
        COUNT += 1
        print('the', COUNT, 'images')

imagepath = 'Y:/2.postgrad/16-liyuanyuan/C_elegent/data_deeplearning/data/train/2.trans_train/'
savepath = 'Y:/2.postgrad/16-liyuanyuan/C_elegent/ARS/seg/data/patch/train/label/'
savepath_r = 'Y:/2.postgrad/16-liyuanyuan/C_elegent/ARS/seg/data/patch/train_resize/label/'
data_form="image"
get_patch_D(imagepath,savepath,savepath_r,data_form)
