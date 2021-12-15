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
from tensorflow.keras.layers import Flatten
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Input
from tensorflow.keras.utils import plot_model


# ToDo: Generate a proper docstring
# This code generates the model that will be loaded by JANE.


class Cerebro:
    @staticmethod
    def conv_activation_pool(layer_input: int, nf: int = 32, sf: int = 3,
                             st: int = 1, act: str = "relu", ps: int = 3,
                             do: float = 0.25, pdd: str = "same", explicit: bool = False):
        if explicit:
            x = Conv1D(filters=nf, kernel_size=sf, strides=st, padding=pdd)(layer_input)
            x = Activation(activation=act)(x)
            x = BatchNormalization()(x)
            x = MaxPooling1D(pool_size=ps)(x)
            x = Dropout(do)(x)
        else:
            x = Conv1D(filters=nf, kernel_size=sf, activation=act, strides=st, padding=pdd)(layer_input)
            x = BatchNormalization()(x)
            x = MaxPooling1D(pool_size=ps)(x)
            x = Dropout(do)(x)
        return x

    @staticmethod
    def dense_act_batchnorm_dropout(inputs: int, neurons: int, act: str = "relu",
                                    dropout: float = 0.5, explicit: bool = False):
        if explicit:
            x = Dense(neurons)(inputs)
            x = Activation(activation=act)(x)
            x = BatchNormalization()(x)
            x = Dropout(dropout)(x)
        else:
            x = Dense(neurons, activation=act)(inputs)
            x = BatchNormalization()(x)
            x = Dropout(dropout)(x)
        return x

    @staticmethod
    def build_input_spec_branch(inputs, number_neurons: int, final_act: str = "relu", explicit: bool = False):
        # CONV => RELU => POOL (I)
        x = Cerebro.conv_activation_pool(inputs, 64, 10, 1, "relu", 3, 0.25, "same", explicit)

        # CONV => RELU => POOL (II)
        x = Cerebro.conv_activation_pool(x, 32, 5, 1, "relu", 3, 0.25, "same", explicit)

        # CONV => RELU => POOL (III)
        x = Cerebro.conv_activation_pool(x, 32, 3, 1, "relu", 3, 0.10, "same", explicit)

        # Output from input spectrum branch
        x = Flatten()(x)
        if explicit:
            x = Dense(number_neurons)(x)
            x = Activation(activation=final_act, name="spectra_intermediate")(x)
        else:
            x = Dense(number_neurons, activation=final_act, name="spectra_intermediate")(x)
        return x

    @staticmethod
    def build_input_magn_branch(inputs: int, number_neurons: int, final_act: str = "relu", explicit: bool = False):
        # 256 neurons, relu, 50% dropout
        x = Cerebro.dense_act_batchnorm_dropout(inputs, 128, "relu", 0.5, explicit)

        # 126 neurons, relu, 25% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 64, "relu", 0.25, explicit)

        # 32 neurons, relu, 10% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 32, "relu", 0.25, explicit)

        # Output from input magnitudes branch
        x = Flatten()(x)
        if explicit:
            x = Dense(number_neurons)(x)
            x = Activation(activation=final_act, name="magnitude_intermediate")(x)
        else:
            x = Dense(number_neurons, activation=final_act, name="magnitude_intermediate")(x)
        return x

    @staticmethod
    def build_output_sfh_branch(inputs: int, number_neurons: int, final_act: str = "linear", explicit: bool = False):
        # 512 neurons, relu, 50% dropout
        x = Cerebro.dense_act_batchnorm_dropout(inputs, 512, "relu", 0.5, explicit)

        # 256 neurons, relu, 25% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 256, "relu", 0.25, explicit)

        # 256 neurons, relu, 25% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 256, "relu", 0.25, explicit)

        # 128 neurons, relu, 10% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 128, "relu", 0.10, explicit)

        # Output from output sfh branch
        if explicit:
            x = Dense(number_neurons)(x)
            x = Activation(activation=final_act, name="sfh_output")(x)
        else:
            x = Dense(number_neurons, activation=final_act, name="sfh_output")(x)
        return x

    @staticmethod
    def build_output_metallicity_branch(inputs: int, number_neurons: int,
                                        final_act: str = "linear", explicit: bool = False):
        # 512 neurons, relu, 50% dropout
        x = Cerebro.dense_act_batchnorm_dropout(inputs, 512, "relu", 0.5, explicit)

        # 256 neurons, relu, 25% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 256, "relu", 0.25, explicit)

        # 256 neurons, relu, 25% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 256, "relu", 0.25, explicit)

        # 128 neurons, relu, 10% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 128, "relu", 0.10, explicit)

        # Output from output sfh branch
        if explicit:
            x = Dense(number_neurons)(x)
            x = Activation(activation=final_act, name="metallicity_output")(x)
        else:
            x = Dense(number_neurons, activation=final_act, name="metallicity_output")(x)
        return x

    @staticmethod
    def build_model(spectra_data_shape: int, magnitudes_data_shape: int,
                    number_neurons_spec: int, number_neurons_magn: int,
                    number_output_sfh: int, number_output_metal: int = None,
                    intermediate_activation: str = "relu", final_activation: str = "linear"):
        if number_output_metal is None:
            number_output_metal = number_output_sfh

        # Input Layers
        input_spec = Input(shape=(spectra_data_shape, 1), name="spectra_input")
        input_magn = Input(shape=(magnitudes_data_shape, 1), name="magnitude_input")

        # Input Branches
        input_spec_branch = Cerebro.build_input_spec_branch(input_spec, number_neurons_spec, intermediate_activation)
        input_magn_branch = Cerebro.build_input_magn_branch(input_magn, number_neurons_magn, intermediate_activation)

        print(input_spec_branch.shape)
        print(input_magn_branch.shape)

        # Concatenate both input branches
        intermediate_concatted = Concatenate(axis=-1)([input_spec_branch, input_magn_branch])

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

    @staticmethod
    def plot(model, filename="testimage.png"):
        plot_model(model, to_file=filename, show_shapes=True, show_layer_names=True)