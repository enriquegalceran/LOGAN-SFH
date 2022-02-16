# -*- coding: utf-8 -*-
# XAVIER
# eXtragalactic Artificial-neural-network Visualizer and identifIER

from keras.models import Model
from keras.layers import BatchNormalization, Conv1D, MaxPooling1D, Activation, \
    Dropout, Concatenate, Flatten, Dense, Input, Lambda
from tensorflow.keras.utils import plot_model
from tensorflow.keras.optimizers import Adam
import keras.backend as kerasbackend
import tensorflow as tf
from typing import Union, Type
from tensorflow.keras import initializers
from RAVEN import print_build_model_parameters


# ToDo: Generate a proper docstring
# This code generates the model that will be loaded by JANE.


class Cerebro:
    @staticmethod
    def conv_activation_pool(layer_input: int, n_filters: int = 32, filter_size: int = 3,
                             stride: int = 1, act: str = "relu", pool_size: int = 3,
                             dropout: float = 0.25, padding: str = "same", explicit: bool = False,
                             kernel_initializer: Union[str, Type[initializers.GlorotUniform]] = "glorot_uniform"):
        if explicit:
            x = Conv1D(filters=n_filters, kernel_size=filter_size,
                       strides=stride, padding=padding,
                       kernel_initializer=kernel_initializer)(layer_input)
            x = Activation(activation=act)(x)
            x = BatchNormalization()(x)
            x = MaxPooling1D(pool_size=pool_size)(x)
            x = Dropout(dropout)(x)
        else:
            x = Conv1D(filters=n_filters, kernel_size=filter_size, activation=act,
                       strides=stride, padding=padding,
                       kernel_initializer=kernel_initializer)(layer_input)
            x = BatchNormalization()(x)
            x = MaxPooling1D(pool_size=pool_size)(x)
            x = Dropout(dropout)(x)
        return x

    @staticmethod
    def dense_act_batchnorm_dropout(inputs: int, neurons: int, act: str = "relu",
                                    dropout: float = 0.5, explicit: bool = False,
                                    kernel_initializer: Union[str, Type[initializers.GlorotUniform]] =
                                    "glorot_uniform"):
        if explicit:
            x = Dense(neurons, kernel_initializer=kernel_initializer)(inputs)
            x = Activation(activation=act)(x)
            x = BatchNormalization()(x)
            x = Dropout(dropout)(x)
        else:
            x = Dense(neurons, kernel_initializer=kernel_initializer, activation=act)(inputs)
            x = BatchNormalization()(x)
            x = Dropout(dropout)(x)
        return x

    @staticmethod
    def build_iterative_branch(inputs, branch_type: str, layers: list = None, number_layers: int = None,
                               neurons_first_layer: int = 128,
                               progression=0.5, filter_size=30,
                               stride=1, act="relu", pool_size=3, dropout=0.25,
                               neurons_last_layer: int = None, explicit: bool = False,
                               output: int = 1, output_neurons: int = 20, final_act: str = "relu",
                               final_layer_name: str = None,
                               kernel_initializer: Union[str, Type[initializers.GlorotUniform]] = "glorot_uniform",
                               ):
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
        :param neurons_first_layer: this parameter represents the number of neurons on a given layer in a dense layer,
                                    or the number of filters in a convolutional layer
        :param progression:
        :param filter_size:
        :param stride:
        :param act:
        :param pool_size:
        :param dropout:
        :param neurons_last_layer:
        :param explicit:
        :param output:              output == 0 -> no output
                                    output == 1 -> output with Flatten (input branches)
                                    output == 2 -> output without Flatten (output branches)
        :param output_neurons:
        :param final_act:
        :param final_layer_name:
        :param kernel_initializer:
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

        filter_size = verify_list(filter_size, "filter_size", number_layers)
        stride = verify_list(stride, "stride", number_layers)
        act = verify_list(act, "act", number_layers)
        pool_size = verify_list(pool_size, "pool_size", number_layers)
        dropout = verify_list(dropout, "dropout", number_layers)
        kernel_initializer = verify_list(kernel_initializer, "kernel_initializer", number_layers)

        # Generate branch layers
        if branch_type.lower() == "cnn":
            for i_layer in range(number_layers):
                if i_layer == 0:
                    x = Cerebro.conv_activation_pool(layer_input=inputs, n_filters=layers[i_layer],
                                                     filter_size=filter_size[i_layer], stride=stride[i_layer],
                                                     act=act[i_layer], pool_size=pool_size[i_layer],
                                                     dropout=dropout[i_layer], padding="same", explicit=explicit,
                                                     kernel_initializer=kernel_initializer[i_layer])
                else:
                    x = Cerebro.conv_activation_pool(layer_input=x, n_filters=layers[i_layer],
                                                     filter_size=filter_size[i_layer], stride=stride[i_layer],
                                                     act=act[i_layer], pool_size=pool_size[i_layer],
                                                     dropout=dropout[i_layer], padding="same", explicit=explicit,
                                                     kernel_initializer=kernel_initializer[i_layer])
        elif branch_type.lower() == "dense":
            for i_layer in range(number_layers):
                if i_layer == 0:
                    x = Cerebro.dense_act_batchnorm_dropout(inputs=inputs, neurons=layers[i_layer], act=act[i_layer],
                                                            dropout=dropout[i_layer], explicit=explicit,
                                                            kernel_initializer=kernel_initializer[i_layer])
                else:
                    x = Cerebro.dense_act_batchnorm_dropout(inputs=x, neurons=layers[i_layer], act=act[i_layer],
                                                            dropout=dropout[i_layer], explicit=explicit,
                                                            kernel_initializer=kernel_initializer[i_layer])
        else:
            raise ValueError("branch_type needs to be 'cnn' or 'dense'.")

        # Output
        if output == 0:
            return x
        elif output == 1:
            x = Flatten()(x)

        if explicit:
            x = Dense(output_neurons)(x)
            x = Activation(activation=final_act, name=final_layer_name)(x)
        else:
            x = Dense(output_neurons, activation=final_act, name=final_layer_name)(x)

        return x

    @staticmethod
    def splitkwargs(spectra_prefix="spect_", magnitude_prefix="magn_", sfh_prefix="sfh_", metal_prefix="metal_",
                    **kwargs):
        """
        Splits kwargs given to build_model depending on the prefix.
        Returns variables with which to update the arguments for each branch.

        :param spectra_prefix:
        :param magnitude_prefix:
        :param sfh_prefix:
        :param metal_prefix:
        :param kwargs:
        :return:
        """

        def separate_by_prefix(args, prefix):
            dictionary = dict()
            for key, value in args.items():
                if prefix in key:
                    dictionary[key[len(prefix):]] = value
            return dictionary

        # For each prefix, filter out the
        spectr_arguments = separate_by_prefix(kwargs, spectra_prefix)
        magnit_arguments = separate_by_prefix(kwargs, magnitude_prefix)
        sfh_arguments = separate_by_prefix(kwargs, sfh_prefix)
        metal_arguments = separate_by_prefix(kwargs, metal_prefix)

        # Remove 'inputs' from dictionary if it was included. This would clash, as it gives multiple instances of
        # the same argument
        spectr_arguments.pop("inputs", None)
        magnit_arguments.pop("inputs", None)
        sfh_arguments.pop("inputs", None)
        metal_arguments.pop("inputs", None)

        return spectr_arguments, magnit_arguments, sfh_arguments, metal_arguments

    @staticmethod
    def build_model(epochs: int = 50, input_mode: str = "single", loss_function_used: str = "SMAPE",
                    loss_function_used_metal: str = None, init_lr: float = 1e-3, loss_weights: tuple = (1.0, 0.8),
                    explicit: bool = False, spectra_data_shape: tuple = (3761, 1),
                    magnitudes_data_shape: tuple = (5, 1), single_data_shape: tuple = (3766, 1),
                    agevector_data_shape=17, verbose: int = 1, **kwargs):

        # Define the default arguments for each branch
        spectr_arguments = {"branch_type": "cnn", "number_layers": 3, "neurons_first_layer": 128,
                            "progression": 0.5, "filter_size": [30, 10, 3], "stride": 1, "act": "relu",
                            "pool_size": 3, "dropout": 0.15, "explicit": explicit,
                            "output": 1, "output_neurons": 128, "final_act": "relu",
                            "final_layer_name": "spectra_intermediate",
                            "kernel_initializer": "glorot_uniform",
                            }
        magn_arguments = {"branch_type": "dense", "number_layers": 2, "neurons_first_layer": 64,
                          "progression": 0.5, "act": "relu",
                          "dropout": 0.25, "explicit": explicit,
                          "output": 1, "output_neurons": 32, "final_act": "relu",
                          "final_layer_name": "magnitude_intermediate",
                          "kernel_initializer": "glorot_uniform",
                          }
        sfh_arguments = {"branch_type": "dense", "layers": [512, 256, 256, 128], "act": "relu",
                         "dropout": 0.15, "explicit": explicit,
                         "output": 2, "output_neurons": agevector_data_shape, "final_act": "relu",
                         "final_layer_name": "sfh_output"}
        metal_arguments = {"branch_type": "dense", "layers": [512, 256, 256, 128], "act": "relu",
                           "dropout": 0.15, "explicit": explicit,
                           "output": 2, "output_neurons": agevector_data_shape, "final_act": "relu",
                           "final_layer_name": "metallicity_output",
                           "kernel_initializer": "glorot_uniform",
                           }

        # Input Layers
        if input_mode == "double":
            input_spec = Input(shape=spectra_data_shape, name="spectra_input")
            input_magn = Input(shape=magnitudes_data_shape, name="magnitude_input")
            input_layer = None  # Defined for compatibility
        elif input_mode == "single":
            input_layer = Input(shape=single_data_shape, name="single_input")
            input_spec = Lambda(lambda x: tf.expand_dims(x[:, :spectra_data_shape[0], 0], -1))(input_layer)
            input_magn = Lambda(lambda x: tf.expand_dims(x[:, spectra_data_shape[0]:, 0], -1))(input_layer)
        else:
            raise ValueError(f"'in_layer' is {input_mode}, and should be either 'single' or 'double'.")

        # Split kwargs for each branch and update
        spectr_kwarg, magnit_kwarg, sfh_kwarg, metal_kwarg = Cerebro.splitkwargs(**kwargs)
        spectr_arguments.update(spectr_kwarg)
        magn_arguments.update(magnit_kwarg)
        sfh_arguments.update(sfh_kwarg)
        metal_arguments.update(metal_kwarg)

        # Input Branches
        input_spec_branch = Cerebro.build_iterative_branch(input_spec, **spectr_arguments)
        input_magn_branch = Cerebro.build_iterative_branch(input_magn, **magn_arguments)

        # Concatenate Both Input Branches
        intermediate_concatted = Concatenate(axis=-1)([input_spec_branch, input_magn_branch])

        # Output Branches
        metal_branch = Cerebro.build_iterative_branch(intermediate_concatted, **metal_arguments)
        sfh_branch = Cerebro.build_iterative_branch(intermediate_concatted, **sfh_arguments)

        # Concatenate Output for single i/o and Define Model
        if input_mode == "double":
            model = Model(
                inputs=[input_spec, input_magn],
                outputs=[metal_branch, sfh_branch],
                name="cerebro_double"
            )
        elif input_mode == "single":
            # Concatenate output
            output_layer = Concatenate(axis=-1)([sfh_branch, metal_branch])
            model = Model(
                inputs=input_layer,
                outputs=output_layer,
                name="cerebro_single"
            )
        else:
            raise ValueError(f"'in_layer' is {input_mode}, and should be either 'single' or 'double'.")

        if loss_function_used == "SMAPE":
            # SMAPE
            # Note: indicate the function, do not call it. (i.e.: 'Cerebro.smape_loss'; NOT 'Crebro.smape_loss()')
            if verbose >= 1:
                print("[INFO] Custom SMAPE is being used.")
            loss_function_used = Cerebro.smape_loss

        if loss_function_used_metal is None:
            loss_function_used_metal = loss_function_used

        if input_mode == "double":
            losses = {
                "sfh_output": loss_function_used,
                "metallicity_output": loss_function_used_metal
            }

            # Loss weight for the two different branches
            loss_weights = {
                "sfh_output": loss_weights[0],
                "metallicity_output": loss_weights[1]
            }
        elif input_mode == "single":
            losses = loss_function_used
            loss_weights = None
        else:
            raise ValueError(f"'in_layer' is {input_mode}, and should be either 'single' or 'double'.")

        # Initialize optimizer and compile the model
        if verbose >= 1:
            print("[INFO] Compiling model...")
        opt = Adam(learning_rate=init_lr, decay=init_lr / epochs)
        model.compile(optimizer=opt, loss=losses, loss_weights=loss_weights, metrics=["accuracy"])

        if verbose == 2:
            print("[INFO] Parameters used to build the model:")
            print_build_model_parameters(spectr_arguments, magn_arguments, sfh_arguments, metal_arguments)

        if verbose >= 3:
            print("[INFO] Model Summary:")
            model.summary()


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
