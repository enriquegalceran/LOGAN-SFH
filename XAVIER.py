# eXtragalactic Artificial-neural-network Visualizer and identifIER

from tensorflow.keras.models import Model
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import Conv1D
from tensorflow.keras.layers import MaxPooling1D
from tensorflow.keras.layers import Activation
from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import Lambda
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Flatten
from tensorflow.keras.layers import Input
import tensorflow as tf

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ToDo: Generate a proper docstring
# This code generates the model that will be loaded by JANE.


class Cerebro:
    @staticmethod
    def conv_activation_pool(layer_input, nf=32, sf=3,
                             st=1, act="relu", ps=3, do=0.25, pd="same"):
        x = Conv1D(filters=nf, kernel_size=sf, strides=st, padding=pd)(layer_input)
        x = Activation(act)(x)
        x = BatchNormalization(axis=1)(x)
        x = MaxPooling1D(pool_size=ps)(x)
        x = Dropout(do)(x)
        return x

    @staticmethod
    def build_input_spectrum_branch(inputs):

        # CONV => RELU => POOL (I)
        x = Cerebro.conv_activation_pool(inputs, 32, 3, 1, "relu", 3, 0.25, "same")

        # CONV => RELU => POOL (II)
        x = Cerebro.conv_activation_pool(x, 32, 3, 1, "relu", 3, 0.25, "same")

        # CONV => RELU => POOL (III)
        x = Cerebro.conv_activation_pool(x, 32, 3, 1, "relu", 3, 0.25, "same")

        return x

    @staticmethod
    def build_input_magnitudes_branch(inputs,):
        pass

    @staticmethod
    def build_output_sfh_branch(inputs,):
        pass

    @staticmethod
    def build_output_metallicity_branch(inputs,):
        pass

    @staticmethod
    def build_model(inputs,):
        pass






