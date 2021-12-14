# -*- coding: utf-8 -*-
# XAVIER
# eXtragalactic Artificial-neural-network Visualizer and identifIER

from tensorflow.keras.models import Model
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import Conv1D
from tensorflow.keras.layers import MaxPooling1D
from tensorflow.keras.layers import Activation
from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import Concatenate
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Input


# ToDo: Generate a proper docstring
# This code generates the model that will be loaded by JANE.


class Cerebro:
    @staticmethod
    def conv_activation_pool(layer_input: int, nf: int = 32, sf: int = 3,
                             st: int = 1, act: str = "relu", ps: int = 3,
                             do: float = 0.25, pdd: str = "same"):
        x = Conv1D(filters=nf, kernel_size=sf, strides=st, padding=pdd)(layer_input)
        x = Activation(act)(x)
        x = BatchNormalization()(x)
        x = MaxPooling1D(pool_size=ps)(x)
        x = Dropout(do)(x)
        return x

    @staticmethod
    def dense_act_batchnorm_dropout(inputs: int, neurons: int, activation: str = "relu", dropout: float = 0.5):
        x = Dense(neurons)(inputs)
        x = Activation(activation)(x)
        x = BatchNormalization()(x)
        x = Dropout(dropout)(x)
        return x

    @staticmethod
    def build_input_spec_branch(inputs: int, number_neurons: int, final_act: str = "relu"):
        # CONV => RELU => POOL (I)
        x = Cerebro.conv_activation_pool(inputs, 32, 3, 1, "relu", 3, 0.25, "same")

        # CONV => RELU => POOL (II)
        x = Cerebro.conv_activation_pool(x, 32, 3, 1, "relu", 3, 0.25, "same")

        # CONV => RELU => POOL (III)
        x = Cerebro.conv_activation_pool(x, 32, 3, 1, "relu", 3, 0.25, "same")

        # Output from input spectrum branch
        x = Dense(number_neurons)(x)
        x = Activation(final_act, name="spectra_intermediate")(x)
        return x

    @staticmethod
    def build_input_magn_branch(inputs: int, number_neurons: int, final_act: str = "relu"):
        # 256 neurons, relu, 50% dropout
        x = Cerebro.dense_act_batchnorm_dropout(inputs, 256, "relu", 0.5)

        # 126 neurons, relu, 25% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 126, "relu", 0.25)

        # Output from input magnitudes branch
        x = Dense(number_neurons)(x)
        x = Activation(final_act, name="magnitude_intermediate")(x)
        return x

    @staticmethod
    def build_output_sfh_branch(inputs: int, number_neurons: int, final_act: str = "linear"):
        # 512 neurons, relu, 50% dropout
        x = Cerebro.dense_act_batchnorm_dropout(inputs, 512, "relu", 0.5)

        # 256 neurons, relu, 25% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 256, "relu", 0.25)

        # 128 neurons, relu, 10% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 128, "relu", 0.10)

        # Output from output sfh branch
        x = Dense(number_neurons)(x)
        x = Activation(final_act, name="sfh_output")(x)
        return x

    @staticmethod
    def build_output_metallicity_branch(inputs: int, number_neurons: int, final_act: str = "linear"):
        # 512 neurons, relu, 50% dropout
        x = Cerebro.dense_act_batchnorm_dropout(inputs, 512, "relu", 0.5)

        # 256 neurons, relu, 25% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 256, "relu", 0.25)

        # 256 neurons, relu, 25% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 256, "relu", 0.25)

        # 128 neurons, relu, 10% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 128, "relu", 0.10)

        # Output from output sfh branch
        x = Dense(number_neurons)(x)
        x = Activation(final_act, name="metallicity_output")(x)
        return x

    @staticmethod
    def build_model(spectra_data_shape: int, magnitudes_data_shape: int,
                    number_neurons_spec: int, number_neurons_magn: int,
                    number_output_metal: int, number_output_sfh: int = None,
                    intermediate_activation: str = "relu", final_activation: str = "linear"):
        if number_output_sfh is None:
            number_output_sfh = number_output_metal

        # Input Layers
        input_spec = Input(shape=spectra_data_shape, name="spectra_input")
        input_magn = Input(shape=magnitudes_data_shape, name="magnitude_input")

        # Input Branches
        input_spec_branch = Cerebro.build_input_spec_branch(input_spec, number_neurons_spec, intermediate_activation)
        input_magn_branch = Cerebro.build_input_magn_branch(input_magn, number_neurons_magn, intermediate_activation)

        # Concatenate both input branches
        intermediate_concatted = Concatenate(axis=0)([input_spec_branch, input_magn_branch])

        metal_branch = Cerebro.build_output_metallicity_branch(intermediate_concatted,
                                                               number_output_metal, final_activation)
        sfh_branch = Cerebro.build_output_sfh_branch(intermediate_concatted,
                                                     number_output_sfh, final_activation)

        model = Model(
            inputs=[input_spec, input_magn],
            outputs=[metal_branch, sfh_branch],
            name="cerebro"
        )

        # ToDo: WIP
        return model
