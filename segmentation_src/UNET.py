import torch
import torch.nn as nn
from torch.autograd import Variable
import torch.nn.functional as F
class Doubleconv(nn.Module):
    def __init__(self,in_ch,out_ch):
        super(Doubleconv, self).__init__()
        self.conv = nn.Sequential(
            nn.Conv3d(in_ch,out_ch,3,padding=1),
            nn.GroupNorm(int(out_ch/2),out_ch),
            #nn.BatchNorm3d(out_ch),
            nn.ReLU(inplace=True),

            nn.Conv3d(out_ch, out_ch, 3, padding=1),
            nn.GroupNorm(int(out_ch/2),out_ch),
            #nn.BatchNorm3d(out_ch),
            nn.ReLU(inplace=True)
        )
    def forward(self, x):
        return  self.conv(x)

class Cas_unet3d(nn.Module):
    def __init__(self,in_ch,out_ch):
        super(Cas_unet3d, self).__init__()
        #encode_layer1  out_ch
        self.conv1 = Doubleconv(in_ch,32)
        self.pool1 = nn.MaxPool3d(2,stride=2,ceil_mode=True)

        #encode_layer2
        self.conv2 =Doubleconv(32,64)
        self.pool2 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        #encode_layer3
        self.conv3 =Doubleconv(64,128)
        self.pool3 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        #encode_layer4
        self.conv4 =Doubleconv(128,256)
        self.pool4 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        # encode_layer5
        self.conv5 = Doubleconv(256, 512)
        self.pool5 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        #decode_upsampling
        self.up6 = nn.ConvTranspose3d(512,256,2,stride=2)
        self.conv6 = Doubleconv(512,256)

        self.up7 = nn.ConvTranspose3d(256, 128, 2, stride=2)
        self.conv7 = Doubleconv(256, 128)

        self.up8 = nn.ConvTranspose3d(128, 64, 2, stride=2)
        self.conv8 = Doubleconv(128, 64)

        self.up9 = nn.ConvTranspose3d(64, 32, 2, stride=2)
        self.conv9 = Doubleconv(64, 32)

        self.conv10 = nn.Conv3d(32,out_ch,1)
        ###############################################

        self.conv2_1 = Doubleconv(out_ch,32)
        self.pool2_1 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        self.conv2_2 = Doubleconv(32, 64)
        self.pool2_2 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        self.conv2_3 = Doubleconv(64, 128)
        self.pool2_3 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        self.conv2_4 = Doubleconv(128, 256)
        self.pool2_4 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        self.conv2_5 = Doubleconv(256, 512)


        self.up2_6 = nn.ConvTranspose3d(512, 256, 2, stride=2)
        self.conv2_6 = Doubleconv(512, 256)

        self.up2_7 = nn.ConvTranspose3d(256, 128, 2, stride=2)
        self.conv2_7 = Doubleconv(256, 128)

        self.up2_8 = nn.ConvTranspose3d(128, 64, 2, stride=2)
        self.conv2_8 = Doubleconv(128, 64)

        self.up2_9 = nn.ConvTranspose3d(64, 32, 2, stride=2)
        self.conv2_9 = Doubleconv(64, 32)

        self.conv2_10 = nn.Conv3d(32, out_ch, 1)

    def forward(self, x):
        c1 = self.conv1(x)
        p1 = self.pool1(c1)

        c2 = self.conv2(p1)
        p2 = self.pool2(c2)

        c3 = self.conv3(p2)
        p3 = self.pool3(c3)

        c4 = self.conv4(p3)
        p4 = self.pool4(c4)

        c5 = self.conv5(p4)

        up_6 =self.up6(c5)
        merge6 = torch.cat([up_6,c4],dim=1)
        c6 = self.conv6(merge6)

        up_7 = self.up7(c6)
        merge7 =  torch.cat([up_7,c3],dim=1)
        c7 = self.conv7(merge7)

        up_8 = self.up8(c7)
        merge8 = torch.cat([up_8, c2], dim=1)
        c8 = self.conv8(merge8)


        up_9 = self.up9(c8)
        merge9 = torch.cat([up_9, c1], dim=1)
        c9 = self.conv9(merge9)

        c10 = self.conv10(c9)

        c2_1 = self.conv2_1(c10)
        p2_1 = self.pool2_1(c2_1)

        c2_2 = self.conv2_2(p2_1)
        p2_2 = self.pool2_2(c2_2)

        c2_3 = self.conv2_3(p2_2)
        p2_3 = self.pool2_3(c2_3)

        c2_4 = self.conv2_4(p2_3)
        p2_4 = self.pool2_4(c2_4)

        c2_5 = self.conv2_5(p2_4)

        up_2_6 = self.up2_6(c5)
        merge2_6 = torch.cat([up_2_6, c2_4], dim=1)
        c2_6 = self.conv2_6(merge2_6)

        up_2_7 = self.up2_7(c2_6)
        merge2_7 = torch.cat([up_2_7, c2_3], dim=1)
        c2_7 = self.conv2_7(merge2_7)

        up_2_8 = self.up2_8(c2_7)
        merge2_8 = torch.cat([up_2_8, c2_2], dim=1)
        c2_8 = self.conv2_8(merge2_8)

        up_2_9 = self.up2_9(c2_8)
        merge2_9 = torch.cat([up_2_9, c2_1], dim=1)
        c2_9 = self.conv2_9(merge2_9)

        c2_10 = self.conv2_10(c2_9)

        return c2_10


class unet3d(nn.Module):
    def __init__(self,in_ch,out_ch):
        super(unet3d, self).__init__()
        #encode_layer1
        self.conv1 =Doubleconv(in_ch,32)
        self.pool1 = nn.MaxPool3d(2,stride=2,ceil_mode=True)

        #encode_layer2
        self.conv2 =Doubleconv(32,64)
        self.pool2 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        #encode_layer3
        self.conv3 =Doubleconv(64,128)
        self.pool3 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        #encode_layer4
        self.conv4 =Doubleconv(128,256)
        self.pool4 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        # encode_layer5
        self.conv5 = Doubleconv(256, 512)
        self.pool5 = nn.MaxPool3d(2, stride=2, ceil_mode=True)

        #decode_upsampling
        self.up6 = nn.ConvTranspose3d(512,256,2,stride=2)
        self.conv6 = Doubleconv(512,256)

        self.up7 = nn.ConvTranspose3d(256, 128, 2, stride=2)
        self.conv7 = Doubleconv(256, 128)

        self.up8 = nn.ConvTranspose3d(128, 64, 2, stride=2)
        self.conv8 = Doubleconv(128, 64)

        self.up9 = nn.ConvTranspose3d(64, 32, 2, stride=2)
        self.conv9 = Doubleconv(64, 32)

        self.conv10 = nn.Conv3d(32,out_ch,1)

    def forward(self, x):
        c1 = self.conv1(x)
        p1 = self.pool1(c1)

        c2 = self.conv2(p1)
        p2 = self.pool2(c2)

        c3 = self.conv3(p2)
        p3 = self.pool3(c3)

        c4 = self.conv4(p3)
        p4 = self.pool4(c4)

        c5 = self.conv5(p4)

        up_6 =self.up6(c5)

        merge6 = torch.cat([up_6,c4],dim=1)
        c6 = self.conv6(merge6)

        up_7 = self.up7(c6)
        merge7 = merge6 = torch.cat([up_7,c3],dim=1)
        c7 = self.conv7(merge7)

        up_8 = self.up8(c7)
        merge8 = torch.cat([up_8, c2], dim=1)
        c8 = self.conv8(merge8)

        up_9 = self.up9(c8)
        merge9 = torch.cat([up_9, c1], dim=1)
        c9 = self.conv9(merge9)

        c10 = self.conv10(c9)

        out = c10
        #out = nn.Sigmoid()(c10)
        return out


# unet = unet3d(3,1)
# print(unet)
