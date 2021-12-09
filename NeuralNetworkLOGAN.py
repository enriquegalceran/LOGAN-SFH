import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# import the necessary packages
from tensorflow.keras.models import Model
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import MaxPooling2D
from tensorflow.keras.layers import Activation
from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import Lambda
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Flatten
from tensorflow.keras.layers import Input
import tensorflow as tf







# This is brought from my TFM
# def CNN():
#     ConvNNet = Sequential()
#     ConvNNet.add(Conv1D(filters = 32, kernel_size = 10,activation='linear',input_shape=(FeaturesTRN.shape[1],1),padding='same', name='Block1CCN'))
#     ConvNNet.add(LeakyReLU(alpha=0.1, name='LeakyReLu1'))
#     ConvNNet.add(MaxPooling1D(pool_size = 4, padding='same', name='MaxPooling1'))
#     ConvNNet.add(Conv1D(filters =64, kernel_size = 5, activation='linear',padding='same', name='Block2CCN'))
#     ConvNNet.add(LeakyReLU(alpha=0.1, name='LeakyReLu2'))
#     ConvNNet.add(MaxPooling1D(pool_size= 2, padding='same', name='MaxPooling2'))
#     ConvNNet.add(Conv1D(filters = 128, kernel_size =3 , activation='linear',padding='same', name='Block3CCN'))
#     ConvNNet.add(LeakyReLU(alpha=0.1, name='LeakyReLu3'))
#     ConvNNet.add(MaxPooling1D(pool_size= 2 ,padding='same', name='MaxPooling3'))
#     ConvNNet.add(Flatten(name='Flatten1'))
#     ConvNNet.add(Dense(128, activation='linear', use_bias = False, name='Dense1'))
#     ConvNNet.add(BatchNormalization(name='BatchNorm1'))
#     ConvNNet.add(LeakyReLU(alpha=0.1, name='LeakyReLu4'))
#     #ConvNNet.add(Dropout(0.1, name='Dropout1'))
#     ConvNNet.add(Dense(num_classes, activation='softmax', name='Dense2'))
#     return ConvNNet





