import time
import numpy as np
import torch
import argparse
import segmentation_src.UNET
from torch import optim
from segmentation_src.dataset import C_Dataset
from torch.utils.data import DataLoader
import os
from torch.optim.lr_scheduler import ReduceLROnPlateau
from segmentation_src.model import dist_loss,dist_loss_mul,color_loss


#os.environ["CUDA_VISIBLE_DEVICES"] = "0"

parser = argparse.ArgumentParser() #创建一个ArgumentParser对象
parser.add_argument('--action', type=str, help='train or test',default = 'train')#添加参数
parser.add_argument('--batch_size', type=int, default=1)
parser.add_argument('--weight', type=str, help='the path of the mode weight file')
parser.add_argument("--lr", type=float,
                    dest="lr", default=1e-4,help="learning rate")
parser.add_argument("--modelpath", type=str,
                    dest="modelpath", default='/media/qulab/720EA7AE0EA76A351/WM/NEW/MODEL/patch_80_c3_7_1.pth',
                    help="frequency of saving models")
args = parser.parse_args()

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#torch.backends.cudnn.enabled = False
lr = args.lr
batch_size = args.batch_size
n_checkpoint = args.checkpoint

model_dir = './MODEL/'
num_epochs = 161
if not os.path.isdir(model_dir):
    os.mkdir(model_dir)
    

def train_model():
    print(torch.cuda.is_available())
    model = UNET.unet3d(1, 3).to(device)
   # model.load_state_dict(torch.load(args.modelpath))
    # 加载数据集
    C_dataset = C_Dataset("./data/seg_sample/train/imgpatch/")
    dataloader = DataLoader(C_dataset, batch_size=batch_size, shuffle=False, num_workers=0)
    
    optimizer = optim.Adam(model.parameters(), lr=lr)
    scheduler = ReduceLROnPlateau(optimizer, 'min',factor=0.8, patience=5)
    min_loss = 100000
    if os.path.exists(model_dir + 'patch_80_c3_9_1_mul.pth'):
        model.load_state_dict(torch.load(args.modelpath))
        print('load old model!')
    
    for epoch in range(num_epochs):
        print('Epoch {}/{}'.format(epoch,num_epochs-1))
        print('-'*10)
        dataset_size = len(dataloader.dataset)
        step = 0
        sum_loss = 0
        for x,y in dataloader:
            start_time = time.time()
            inputs = x.to(device)
            labels = y.to(device)
            y1 = (y*y).sum(axis=1)
            count = np.where(y1>0)
            outputs = model(inputs)
            print('y_pre:',outputs.max())
            print('labels:',labels.max())
            loss1 = 7 *dist_loss_mul(labels,outputs)
            loss2 = color_loss(labels,outputs,y1)
            print("loss1=",loss1)
            print("loss2=",loss2)
            loss = loss1 + loss2
            lo=loss.item()
            sum_loss=lo + sum_loss
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            print('[', "{0:.3f}".format((time.time() - start_time)), 'seconds]', 'After', step + 1,
                  'steps, the training loss is ', loss)
            step += 1
        aver_loss =  sum_loss / dataset_size
        if aver_loss < min_loss:
            min_loss = aver_loss
            print("save model")
            torch.save(model.state_dict(), model_dir + '/' + 'patch_80_c3_9_1_mul' + '.pth')
        print("%d epoch aver_loss:%10f"%(epoch,aver_loss))
        print("%d epoch %10f" % (epoch, optimizer.param_groups[0]['lr']))
        scheduler.step(aver_loss)

    return model
train_model()


