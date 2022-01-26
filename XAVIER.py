# -*- coding: utf-8 -*-
# XAVIER
# eXtragalactic Artificial-neural-network Visualizer and identifIER

from tensorflow.keras.models import Model
from tensorflow.keras.layers import BatchNormalization, Conv1D, MaxPooling1D, Activation, \
    Dropout, Concatenate, Flatten, Dense, Input
from tensorflow.keras.utils import plot_model
import tensorflow.keras.backend as kerasbackend


# ToDo: Generate a proper docstring
# This code generates the model that will be loaded by JANE.


class Cerebro:
    @staticmethod
    def conv_activation_pool(layer_input: int, n_filters: int = 32, filter_size: int = 3,
                             stride: int = 1, act: str = "relu", pool_size: int = 3,
                             dropout: float = 0.25, padding: str = "same", explicit: bool = False):
        if explicit:
            x = Conv1D(filters=n_filters, kernel_size=filter_size,
                       strides=stride, padding=padding)(layer_input)
            x = Activation(activation=act)(x)
            x = BatchNormalization()(x)
            x = MaxPooling1D(pool_size=pool_size)(x)
            x = Dropout(dropout)(x)
        else:
            x = Conv1D(filters=n_filters, kernel_size=filter_size, activation=act,
                       strides=stride, padding=padding)(layer_input)
            x = BatchNormalization()(x)
            x = MaxPooling1D(pool_size=pool_size)(x)
            x = Dropout(dropout)(x)
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
    def build_iterative_branch(inputs, branch_type: str, layers: list = None, number_layers: int = None,
                               neurons_first_layer: int = 128,
                               progression=0.5, filter_size=30,
                               stride=1, act="relu", pool_size=3, dropout=0.25,
                               neurons_last_layer: int = None, explicit: bool = False):
        """
        Generate a branch with lists. Used for CV.
        branch_type needs to be 'cnn' or 'dense'.
        For the variables [filter_size, stride, act, pool_size, dropout], if the input is a list, the list will be used,
        if it is a single value, a list with each value will be used.
        Either layers or number_layers needs to be provided.

        :param inputs:
        :param branch_type:
        :param layers:
        :param number_layers:
        :param neurons_first_layer:
        :param progression:
        :param filter_size:
        :param stride:
        :param act:
        :param pool_size:
        :param dropout:
        :param neurons_last_layer:
        :param explicit:
        :return:
        """
        # Initialize x
        x = None

        # Define number of layers
        if layers is not None:
            number_layers = len(layers)
        elif number_layers is not None:
            layers = [int(neurons_first_layer * (progression ** i)) for i in range(1, number_layers + 1)]
        else:
            raise ValueError("Either layers or number_layers needs to be provided.")
        # Las layer can be overwritten
        if neurons_last_layer is not None:
            layers[-1] = neurons_last_layer

        # Verify that input variables are lists or single elements.
        def verify_list(parameter, parameter_str, n_layers):
            if type(parameter) is list:
                assert len(parameter) == n_layers, f"Length of {parameter_str} is inconsistent"
            else:
                parameter = [parameter for _ in range(n_layers)]
            return parameter
        filter_size = verify_list(filter_size, "kernel_size", number_layers)
        stride = verify_list(stride, "stride", number_layers)
        act = verify_list(act, "act", number_layers)
        pool_size = verify_list(pool_size, "pool_size", number_layers)
        dropout = verify_list(dropout, "dropout", number_layers)

        # Generate branch layers
        if branch_type.lower() == "cnn":
            for i_layer in range(number_layers):
                if i_layer == 0:
                    x = Cerebro.conv_activation_pool(layer_input=inputs, n_filters=layers[i_layer],
                                                     filter_size=filter_size[i_layer], stride=stride[i_layer],
                                                     act=act[i_layer], pool_size=pool_size[i_layer],
                                                     dropout=dropout[i_layer], padding="same", explicit=explicit)
                else:
                    x = Cerebro.conv_activation_pool(layer_input=x, n_filters=layers[i_layer],
                                                     filter_size=filter_size[i_layer], stride=stride[i_layer],
                                                     act=act[i_layer], pool_size=pool_size[i_layer],
                                                     dropout=dropout[i_layer], padding="same", explicit=explicit)
        elif branch_type.lower() == "dense":
            for i_layer in range(number_layers):
                if i_layer == 0:
                    x = Cerebro.dense_act_batchnorm_dropout(inputs=inputs, neurons=layers[i_layer], act=act[i_layer],
                                                            dropout=dropout[i_layer], explicit=explicit)
                else:
                    x = Cerebro.dense_act_batchnorm_dropout(inputs=x, neurons=layers[i_layer], act=act[i_layer],
                                                            dropout=dropout[i_layer], explicit=explicit)
        else:
            raise ValueError("branch_type needs to be 'cnn' or 'dense'.")

        return x


    @staticmethod
    def build_input_spec_branch(inputs, number_output_neurons_spect: int,
                                final_act: str = "relu", explicit: bool = False):
        # CONV => RELU => POOL (I)
        x = Cerebro.conv_activation_pool(inputs, 128, 30, 1, "relu", 3, 0.25, "same", explicit)

        # CONV => RELU => POOL (II)
        x = Cerebro.conv_activation_pool(x, 64, 10, 1, "relu", 3, 0.25, "same", explicit)

        # CONV => RELU => POOL (III)
        x = Cerebro.conv_activation_pool(x, 32, 3, 1, "relu", 3, 0.10, "same", explicit)

        # Output from input spectrum branch
        x = Flatten()(x)
        if explicit:
            x = Dense(number_output_neurons_spect)(x)
            x = Activation(activation=final_act, name="spectra_intermediate")(x)
        else:
            x = Dense(number_output_neurons_spect, activation=final_act, name="spectra_intermediate")(x)
        return x

    @staticmethod
    def build_input_magn_branch(inputs: int, number_output_neurons_mag: int,
                                final_act: str = "relu", explicit: bool = False):

        # 32 neurons, relu, 25% dropout
        x = Cerebro.dense_act_batchnorm_dropout(inputs, 32, "relu", 0.25, explicit)

        # 16 neurons, relu, 10% dropout
        x = Cerebro.dense_act_batchnorm_dropout(x, 16, "relu", 0.25, explicit)

        # Output from input magnitudes branch
        x = Flatten()(x)
        if explicit:
            x = Dense(number_output_neurons_mag)(x)
            x = Activation(activation=final_act, name="magnitude_intermediate")(x)
        else:
            x = Dense(number_output_neurons_mag, activation=final_act, name="magnitude_intermediate")(x)
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
                    intermediate_activation: str = "relu", final_activation: str = "linear", explicit: bool = False):
        if number_output_metal is None:
            number_output_metal = number_output_sfh

        # Input Layers
        input_spec = Input(shape=(spectra_data_shape, 1), name="spectra_input")
        input_magn = Input(shape=(magnitudes_data_shape, 1), name="magnitude_input")

        # Input Branches
        input_spec_branch = Cerebro.build_input_spec_branch(input_spec, number_neurons_spec,
                                                            intermediate_activation, explicit)
        input_magn_branch = Cerebro.build_input_magn_branch(input_magn, number_neurons_magn,
                                                            intermediate_activation, explicit)

        # Concatenate both input branches
        intermediate_concatted = Concatenate(axis=-1)([input_spec_branch, input_magn_branch])

        metal_branch = Cerebro.build_output_metallicity_branch(intermediate_concatted,
                                                               number_output_metal, final_activation, explicit)
        sfh_branch = Cerebro.build_output_sfh_branch(intermediate_concatted,
                                                     number_output_sfh, final_activation, explicit)

        model = Model(
            inputs=[input_spec, input_magn],
            outputs=[metal_branch, sfh_branch],
            name="cerebro"
        )

        return model

    @staticmethod
    def graph(model, filename="testimage.png"):
        plot_model(model, to_file=filename, show_shapes=True, show_layer_names=True)

    @staticmethod
    def smape_loss(y_true, y_pred):
        epsilon = 0.1
        summ = kerasbackend.maximum(kerasbackend.abs(y_true) + kerasbackend.abs(y_pred) + epsilon, 0.5 + epsilon)
        smape = kerasbackend.abs(y_pred - y_true) / summ * 2.0
        return smape
