#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import numpy as np
from random import shuffle
import wavio
import wave
import subprocess
import time
import argparse
import random
import pathlib


# In[2]:


sys.path.append("../../src/")
sys.path.append("../../")


# In[3]:


from common.tlopts import display_info;
import common.utils as U;
from Libs.SharedLibs import *


# In[4]:


from Libs.SharedLibs import getFileList;
from Libs.datetime_util import genDataTimeStr;


# In[5]:


# wav_src_dir = "../../datasets/CurrentUse/wav_files/"
wav_root_src_dir = "../../datasets/CurrentUse/wav_files/Multi_Fold/fold2/"
version_string = "version14_multifold_office"
fold_name = "fold1"
model_input_length = 20150;
wav_sampling_rate = 20000;


# ## acdnet original dataset processing codes

# In[7]:


# def convert_sr(src_path, dst_path, sr):
#     print('* {} -> {}'.format(src_path, dst_path))
#     if not os.path.exists(dst_path):
#         os.mkdir(dst_path);
#     for src_file in sorted(glob.glob(os.path.join(src_path, '*.wav'))):
#         dst_file = src_file.replace(src_path, dst_path);
#         subprocess.call('ffmpeg -i {} -ac 1 -ar {} -loglevel error -y {}'.format(
#             src_file, sr, dst_file), shell=True);


# In[8]:


# def preprocess_dataset(to_process_dir):
#     wavlst = getFileList(to_process_dir);
#     save_dir = os.path.join(to_process_dir, '20K');
#     if not os.path.exists(save_dir):
#         os.mkdir(save_dir)
#     sr = 20000;
#     for src_file in wavlst:
#         print(f"original scc_file:{src_file}")
#         wav_name = os.path.basename(src_file)[:-4];
#         wav_name = "{}_20K.wav".format(wav_name);
#         save_path = os.path.join(save_dir, wav_name)
#         print(f"save dir:{save_path}")
#         framerate = 0;
#         with wave.open(src_file, 'rb') as f:
#             framerate = f.getframerate();
#         if framerate != 20000:
#             subprocess.call('ffmpeg -i {} -ac 1 -ar {} -loglevel error -y {}'.format(
#                 src_file, sr, save_path), shell=True);
#             print(f"converted {src_file} sampling rate from {framerate} to 20K")


# ## Generating Single-Fold Training DataSet

# In[9]:


# def create_base_train_dataset(upper_level_dir=None, fold_dirs=None, p_classes=None, n_classes=None, export_path=None):
def create_base_train_dataset(upper_level_dir=None, p_classes=None, n_classes=None, export_path=None):
    total_counter = 0;
    p_counter = 0;
    n_counter = 0;
    train_dataset = {};
    dict_key = 'fold2';
    # for fold in fold_dirs:
    train_dataset[dict_key] = {}
    train_sounds = []
    train_labels = []
    for t in p_classes:
    #Dealing with positive wav files
        p_current_dir = os.path.join(upper_level_dir,'positive', t);
        print(f"work on dir:{p_current_dir}");
        lbl = t[t.rfind('_')+1:]
        tmp_list = getFileList(p_current_dir)
        for f in tmp_list:
            sound = wavio.read(f).data.T[0]
            start = sound.nonzero()[0].min()
            end = sound.nonzero()[0].max()
            sound = sound[start: end + 1]  # Remove silent sections
            train_sounds.append(sound)
            train_labels.append(lbl)
            total_counter += 1;
            p_counter += 1;
            
    # n_classes_root_dir = os.path.join(upper_level_dir,'negative')
    for c in n_classes:
        n_current_dir = os.path.join(upper_level_dir,'negative',c)
        print(f"work on dir:{n_current_dir}");
        n_lbl = 99;#c[:c.find('_')]
        tmp_list2 = getFileList(n_current_dir)
        for f in tmp_list2:
            sound = wavio.read(f).data.T[0]
            start = sound.nonzero()[0].min()
            end = sound.nonzero()[0].max()
            sound = sound[start: end + 1]  # Remove silent sections
            train_sounds.append(sound)
            train_labels.append(n_lbl)
            total_counter += 1;
            n_counter += 1;

    train_dataset[dict_key]['sounds'] = train_sounds
    train_dataset[dict_key]['labels'] = train_labels
    np.savez(export_path, **train_dataset)
    print(f"Training Data is generated and save at {export_path}")
    print(f"Total wav files for trainset is {total_counter}");
    print(f"Total positive wav files for trainset is {p_counter}");
    print(f"Total negative wav files for trainset is {n_counter}");


# In[10]:


# p_classes = ["alarm_52","help_mandrain_71"]
def main():
    upper_level_dir = "../../datasets/CurrentUse/wav_files/Multi_Fold/fold2/train/"
    save_dir = os.path.join("../CurrentUse/generated_datasets/train/{}_{}",version_string,fold_name);
    if not pathlib.Path(save_dir).is_dir():
        os.makedirs(save_dir);
    output_path = os.path.join(save_dir,"{}_train_{}.npz".format(fold_name,genDataTimeStr()));
    p_dirs = ['alarm_52','moaning_56','help_eng_71'];#[f for f in os.listdir(upper_level_dir+"positive")]
    n_dirs = [f for f in os.listdir(upper_level_dir+"negative")]
    print(f"positive classes:{p_dirs}")
    print("negative classes:")
    count1 = 0
    for c in n_dirs:
        count1 += 1
        print(f"{count1}:{c}")
    dataset = create_base_train_dataset(upper_level_dir=upper_level_dir, p_classes=p_dirs, n_classes=n_dirs, export_path=output_path);


# In[11]:


main()


# In[14]:


ds = np.load("../CurrentUse/generated_datasets/train/version12_4class_onesecond_input_office/single_fold_train_20240619110117.npz", allow_pickle=True);
print(ds['fold1'].item()['labels'])


# ## Generating Validation DataSet

# In[ ]:





# In[15]:


class ValGenerator():
    #Generates data for Keras
    def __init__(self, samples, labels, options):
        random.seed(42);
        #Initialization
        print(len(samples));
        self.data = [(samples[i], labels[i]) for i in range (0, len(samples))];
        self.opt = options;
        self.batch_size = len(samples);#88;#options.batchSize // options.nCrops;
        print(f"batch_size:{self.batch_size}");
        self.preprocess_funcs = self.preprocess_setup();
        self.map_dict= {
            '52':1, #alarm
            '56':2, #moaning
            '71':3, #help(english)
            '99':4, #other_sounds 
        };

    def get_data(self):
        #Generate one batch of data
        x, y = self.generate();
        x = np.expand_dims(x, axis=1)
        x = np.expand_dims(x, axis=3)
        # print(x.shape);
        # print(y.shape);
        return x, y

    def generate(self):
        #Generates data containing batch_size samples
        sounds = [];
        labels = [];
        indexes = None;
        for i in range(self.batch_size):
            sound, target = self.data[i];
            target = self.map_dict[str(target)] - 1;
            sound = self.preprocess(sound).astype(np.float32)
            # print(sound)
            label = np.zeros((self.opt.nCrops, self.opt.nClasses));
            label[:,target] = 1;
            print(f"nCrops:{self.opt.nCrops}, nClasses:{self.opt.nClasses}")
            sounds.append(sound);
            labels.append(label);
        """
        #dtype="object" for ValueError: setting an array element with a sequence. 
        The requested array has an inhomogeneous shape after 1 dimensions. 
        The detected shape was (58,) + inhomogeneous part.
        """
        sounds = np.asarray(sounds,dtype="object")
        # expand_sounds = np.expand_dims(np.asarray(sounds,dtype="object"),axis=1); 
        labels = np.asarray(labels);
        # print(f"shape of sounds:{expand_sounds.shape}")
        sounds = sounds.reshape(sounds.shape[0]*sounds.shape[1], sounds.shape[2]);
        labels = labels.reshape(labels.shape[0]*labels.shape[1], labels.shape[2]);
        return sounds, labels;

    def preprocess_setup(self):
        funcs = []
        funcs += [U.padding(self.opt.inputLength // 2),
                  U.normalize(32768.0),
                  U.multi_crop(self.opt.inputLength, 2)] # we use single crop here.

        return funcs

    def preprocess(self, sound):
        for f in self.preprocess_funcs:
            sound = f(sound)

        return sound;


# In[16]:


def getOpts():
    parser = argparse.ArgumentParser(description='Transfer Learning for ACDNet');
    parser.add_argument('--netType', default='ACDNet_TL_Model_Extend',  required=False);
    parser.add_argument('--data', default='../datasets/processed/',  required=False);
    parser.add_argument('--dataset', required=False, default='uec_iot', choices=['10']);
    parser.add_argument('--BC', default=True, action='store_true', help='BC learning');
    parser.add_argument('--strongAugment', default=True,  action='store_true', help='Add scale and gain augmentation');
    #在ipynb中，不能使用parser.parse，要改用parser.parse_known_args()
    opt, unknown = parser.parse_known_args();

    """
    current best setting for accuracy: 96.5
    opt.batchSize = 64;
    opt.LR = 0.1;
    opt.weightDecay = 5e-4;
    opt.momentum = 0.9;
    opt.schedule = [0.3, 0.5, 0.9];
    """
    #Leqarning settings
    opt.batchSize = 64;
    opt.LR = 0.1;
    opt.weightDecay = 5e-4;#1e-2;#5e-3;#5e-4;
    opt.momentum = 0.9;
    opt.nEpochs = 800;
    opt.schedule = [0.3, 0.5, 0.8];
    opt.warmup = 10;
    opt.device = 'cpu';
    # if torch.backends.mps.is_available():
    #     opt.device="mps"; #for apple m2 gpu
    # elif torch.cuda.is_available():
    #     opt.device="cuda:0"; #for nVidia gpu
    # else:
    #     opt.device="cpu"
    print(f"***Use device:{opt.device}");
    # opt.device = torch.device("cuda:0" if  else "cpu");
    #Basic Net Settings
    opt.nClasses = 4#50;
    opt.nFolds = 1;
    opt.splits = [i for i in range(1, opt.nFolds + 1)];
    opt.sr = 20000;
    opt.inputLength = model_input_length;#20150;#30225;
    #Test data
    opt.nCrops = 2;
    opt.TLAcdnetConfig = [8,64,32,64,64,128,128,256,256,512,512,2];
    return opt


# In[17]:


def create_val_dataset_src_npy(dict_key='fold1', upper_level_dir=None, p_classes=None, n_classes=None, export_path=None):
    train_dataset = {};
    # dict_key = 'fold1';
    # for fold in fold_dirs:
    train_dataset[dict_key] = {}
    train_sounds = []
    train_labels = []
    total_counter = 0;
    p_counter = 0;
    n_counter = 0;
    for t in p_classes:
    #Dealing with positive wav files
        p_current_dir = os.path.join(upper_level_dir,'positive', t);
        print(f"work on dir:{p_current_dir}");
        lbl = t[t.rfind('_')+1:]
        tmp_list = getFileList(p_current_dir)
        print(f"class:{t}, item number:{len(tmp_list)}")
        for f in tmp_list:
            sound = wavio.read(f).data.T[0]
            start = sound.nonzero()[0].min()
            end = sound.nonzero()[0].max()
            sound = sound[start: end + 1]  # Remove silent sections
            train_sounds.append(sound)
            train_labels.append(lbl)
            total_counter += 1;
            p_counter += 1;
    for c in n_classes:
        n_current_dir = os.path.join(upper_level_dir,'negative',c)
        print(f"work on dir:{n_current_dir}");
        n_lbl = 99;#c[:c.find('_')]
        tmp_list2 = getFileList(n_current_dir)
        print(f"class:{c}, item number:{len(tmp_list2)}")
        for f in tmp_list2:
            sound = wavio.read(f).data.T[0]
            start = sound.nonzero()[0].min()
            end = sound.nonzero()[0].max()
            sound = sound[start: end + 1]  # Remove silent sections
            train_sounds.append(sound)
            train_labels.append(n_lbl)
            total_counter += 1;
            n_counter += 1;

    train_dataset[dict_key]['sounds'] = train_sounds
    train_dataset[dict_key]['labels'] = train_labels
    np.savez(export_path, **train_dataset)
    print(f"Test Data is generated and save at {export_path}")
    print(f"Total wav files for trainset is {total_counter}");
    print(f"Total positive wav files for trainset is {p_counter}");
    print(f"Total negative wav files for trainset is {n_counter}");


# In[18]:


def val_main1():
    val_upper_level_dir = "../CurrentUse/wav_files/Single_Fold/val/"
    p_dirs = ['alarm_52','moaning_56','help_eng_71'];#[f for f in os.listdir(upper_level_dir+"positive")]
    n_dirs = [f for f in os.listdir(val_upper_level_dir+"negative")]
    save_path = os.path.join("../CurrentUse/generated_datasets/val/",version_string);
    if not pathlib.Path(save_path).is_dir():
        os.makedirs(save_path)
    output_fullname = os.path.join("../CurrentUse/generated_datasets/val/", version_string, "single_fold_val_src_{}.npz".format(genDataTimeStr()));
    print(f"positive classes:{p_dirs}")
    print("negative directories:")
    count1 = 0
    for c in n_dirs:
        count1 += 1
        print(f"{count1}:{c}")
    create_val_dataset_src_npy(upper_level_dir=val_upper_level_dir,p_classes=p_dirs,n_classes=n_dirs,export_path=output_fullname);


# In[19]:


val_main1()


# In[18]:


val_dataset = np.load('../../../uec_iot_ai_models_datasets/',allow_pickle=True);


# In[19]:


# print(len(val_dataset['fold1'].item()['labels']))
# print(val_dataset['fold1'])


# In[21]:


def create_test_compress_npz(val_src_sounds=None, val_src_labels=None, export_path=None):
    opt = getOpts();#opts.parse();
    # display_info(opt);
    # opt.batchSize=88;
    opt.nCrops = 2;
    opt.nClasses= 4;
    opt.sr = 20000;
    opt.inputLength = 20150;#30225;
    val_sounds = [];
    val_labels = [];
   
    start_time = time.perf_counter();
    val_sounds.extend(val_src_sounds);
    val_labels.extend(val_src_labels);
    print(f"len of val_sounds:{len(val_sounds)}, len of val_labels:{len(val_labels)}")
    
    valGen = ValGenerator(val_src_sounds, val_src_labels, opt);
    valX, valY = valGen.get_data();

    np.savez_compressed(export_path, x=valX, y=valY);
    print('compressed npz generated with\n  shape x:{}\n  y:{}\n  took {:.2f} secs'.format(valX.shape, valY.shape, time.perf_counter()-start_time));
    sys.stdout.flush();


# In[ ]:





# In[22]:


src_val_data_npz = "../CurrentUse/generated_datasets/val/version12_4class_onesecond_input_office/single_fold_val_src_20240619111705.npz"
val_data = np.load(src_val_data_npz, allow_pickle=True);
dest_path = "../CurrentUse/generated_datasets/val/{}/final_single_val_{}.npz".format(version_string,genDataTimeStr());
#call to create validation set
print(f"dest_path is {dest_path}")
_sounds = val_data['fold1'].item()['sounds']
_labels = val_data['fold1'].item()['labels']
create_test_compress_npz(_sounds,_labels,dest_path)


# In[23]:


val_data = np.load("../CurrentUse/generated_datasets/val/version12_4class_onesecond_input_office/final_single_val_20240619112602.npz",allow_pickle=True)


# In[24]:


print(val_data['y'])


# In[ ]:





# In[14]:


# test_list = [52 for _ in range(87)]


# In[15]:


# print(test_list)


# In[16]:


# print(len(test_list))


# In[8]:


# test_lbl = 'alarm_52'
# print(test_lbl[test_lbl.rfind('_')+1:])


# In[7]:


# test_lbl2 = '12_rainfall'
# print(test_lbl2[:test_lbl2.find('_')])


# In[ ]:




