# import the necessary packages
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import train_test_split
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


from tensorflow.keras.optimizers import Adam
from tensorflow.keras.preprocessing.image import img_to_array
from sklearn.preprocessing import LabelBinarizer
from sklearn.model_selection import train_test_split
from imutils import paths
import matplotlib.pyplot as plt
import numpy as np
import argparse
import random
import pickle
import cv2
import os


# https://www.pyimagesearch.com/2018/06/04/keras-multiple-outputs-and-multiple-losses/?_ga=2.184316518.1970541148.1638269733-1158531756.1638269733
class FashionNet:
    @staticmethod
    def build_category_branch(inputs, numCategories, finalAct="softmax", chanDim=-1):
        # utilize a lambda layer to convert the 3 channel input to a
        # grayscale representation
        x = Lambda(lambda c: tf.image.rgb_to_grayscale(c))(inputs)

        # CONV => RELU => POOL
        x = Conv2D(32, (3, 3), padding="same")(x)
        x = Activation("relu")(x)
        x = BatchNormalization(axis=chanDim)(x)
        x = MaxPooling2D(pool_size=(3, 3))(x)
        x = Dropout(0.25)(x)

        # (CONV => RELU) * 2 => POOL
        x = Conv2D(64, (3, 3), padding="same")(x)
        x = Activation("relu")(x)
        x = BatchNormalization(axis=chanDim)(x)
        x = Conv2D(64, (3, 3), padding="same")(x)
        x = Activation("relu")(x)
        x = BatchNormalization(axis=chanDim)(x)
        x = MaxPooling2D(pool_size=(2, 2))(x)
        x = Dropout(0.25)(x)

        # (CONV => RELU) * 2 => POOL
        x = Conv2D(128, (3, 3), padding="same")(x)
        x = Activation("relu")(x)
        x = BatchNormalization(axis=chanDim)(x)
        x = Conv2D(128, (3, 3), padding="same")(x)
        x = Activation("relu")(x)
        x = BatchNormalization(axis=chanDim)(x)
        x = MaxPooling2D(pool_size=(2, 2))(x)
        x = Dropout(0.25)(x)

        # define a branch of output layers for the number of different
        # clothing categories (i.e., shirts, jeans, dresses, etc.)
        x = Flatten()(x)
        x = Dense(256)(x)
        x = Activation("relu")(x)
        x = BatchNormalization()(x)
        x = Dropout(0.5)(x)
        x = Dense(numCategories)(x)
        x = Activation(finalAct, name="category_output")(x)

        # return the category prediction sub-network
        return x

    @staticmethod
    def build_color_branch(inputs, numColors, finalAct="softmax",
                           chanDim=-1):
        # CONV => RELU => POOL
        x = Conv2D(16, (3, 3), padding="same")(inputs)
        x = Activation("relu")(x)
        x = BatchNormalization(axis=chanDim)(x)
        x = MaxPooling2D(pool_size=(3, 3))(x)
        x = Dropout(0.25)(x)

        # CONV => RELU => POOL
        x = Conv2D(32, (3, 3), padding="same")(x)
        x = Activation("relu")(x)
        x = BatchNormalization(axis=chanDim)(x)
        x = MaxPooling2D(pool_size=(2, 2))(x)
        x = Dropout(0.25)(x)

        # CONV => RELU => POOL
        x = Conv2D(32, (3, 3), padding="same")(x)
        x = Activation("relu")(x)
        x = BatchNormalization(axis=chanDim)(x)
        x = MaxPooling2D(pool_size=(2, 2))(x)
        x = Dropout(0.25)(x)

        # define a branch of output layers for the number of different
        # colors (i.e., red, black, blue, etc.)
        x = Flatten()(x)
        x = Dense(128)(x)
        x = Activation("relu")(x)
        x = BatchNormalization()(x)
        x = Dropout(0.5)(x)
        x = Dense(numColors)(x)
        x = Activation(finalAct, name="color_output")(x)

        # return the color prediction sub-network
        return x

    @staticmethod
    def build(width, height, numCategories, numColors,
              finalAct="softmax"):
        # initialize the input shape and channel dimension (this code
        # assumes you are using TensorFlow which utilizes channels
        # last ordering)
        inputShape = (height, width, 3)
        chanDim = -1

        # construct both the "category" and "color" sub-networks
        inputs = Input(shape=inputShape)
        categoryBranch = FashionNet.build_category_branch(inputs,
                                                          numCategories, finalAct=finalAct, chanDim=chanDim)
        colorBranch = FashionNet.build_color_branch(inputs,
                                                    numColors, finalAct=finalAct, chanDim=chanDim)

        # create the model using our input (the batch of images) and
        # two separate outputs -- one for the clothing category
        # branch and another for the color branch, respectively
        model = Model(
            inputs=inputs,
            outputs=[categoryBranch, colorBranch],
            name="fashionnet")
        
        # return the constructed network architecture
        return model


# if False:
#     from tensorflow.keras.models import Sequential
#     from numpy import loadtxt
#     import matplotlib.pyplot as plt
#     import visualkeras
#     from tensorflow.keras.layers import Dense, Dropout, Flatten, Activation, Conv1D, MaxPooling1D, BatchNormalization
#     from tensorflow.keras.layers import LeakyReLU
#
#     # load the dataset
#     dataset = loadtxt('NeuralNetworkData/pima-indians-diabetes.csv', delimiter=',')
#     # split into input (X) and output (y) variables
#     X = dataset[:, 0:8]
#     y = dataset[:, 8]
#     # print(X)
#     # print(y)
#
#     # define the keras model
#     model = Sequential()
#     model.add(Dense(12, input_dim=8, activation='relu'))
#     model.add(Dense(8, activation='relu'))
#     model.add(Dense(1, activation='sigmoid'))
#
#     # compile the keras model
#     model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
#
#     print("----------------___________________------------------")
#     model.summary()
#
#     model.save("model.h5")
#     model2 = CNN_tfm()
#     model2.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
#     model2.save("model2.h5")
#
#     # plot_model(model, to_file="model.png")
#     exit()
#
#     # fit the keras model on the dataset
#     output_model = model.fit(X, y, epochs=200, batch_size=10)
#
#     plt.semilogy(output_model.history['loss'], "b-")
#     plt.semilogy(output_model.history['accuracy'], "r-")
#
#     # evaluate the keras model
#     _, accuracy = model.evaluate(X, y)
#     print('Accuracy: %.2f' % (accuracy * 100))
#
#     ###################################
#     ###################################
#     ###################################
#
#     # make probability predictions with the model
#     predictions = model.predict(X)
#     # round predictions
#     rounded = [round(x[0]) for x in predictions]
#     # summarize the first 5 cases
#     for i in range(5):
#         print('%s => %d (expected %d)' % (X[i].tolist(), predictions[i], y[i]))
#
#     plt.show()




def main():
    # initialize the number of epochs to train for, initial learning rate,
    # batch size, and image dimensions
    EPOCHS = 50
    INIT_LR = 1e-3
    BS = 32
    IMAGE_DIMS = (96, 96, 3)

    # initialize our FashionNet multi-output network
    model = FashionNet.build(96, 96,
                             numCategories=5,
                             numColors=3,
                             finalAct="softmax")
    model.summary()
    model.save("testmodel.h5")
    print("END HERE")
    exit()
    #
    # # define two dictionaries: one that specifies the loss method for
    # # each output of the network along with a second dictionary that
    # # specifies the weight per loss
    # losses = {
    #     "category_output": "categorical_crossentropy",
    #     "color_output": "categorical_crossentropy",
    # }
    # # Ponderar la importancia de cada una de las dos categorias!
    # lossWeights = {"category_output": 1.0, "color_output": 1.0}
    #
    # # initialize the optimizer and compile the model
    # print("[INFO] compiling model...")
    # opt = Adam(lr=INIT_LR, decay=INIT_LR / EPOCHS)
    # model.compile(optimizer=opt, loss=losses, loss_weights=lossWeights,
    #               metrics=["accuracy"])
    #
    # # partition the data into training and testing splits using 80% of
    # # the data for training and the remaining 20% for testing
    # # Aqui se le puede meter cualquier numero de funciones. HAY QUE ELIMINAR EL RANDOM_STATE CUANDO FUNCIONE!
    # split = train_test_split(data, categoryLabels, colorLabels,
    #                          test_size=0.2, random_state=42)
    # (trainX, testX, trainCategoryY, testCategoryY,
    #  trainColorY, testColorY) = split
    # # ---------------------------------------------------------------------------------------
    # # train the network to perform multi-output classification
    # H = model.fit(x=trainX,
    #               y={"category_output": trainCategoryY, "color_output": trainColorY},
    #               validation_data=(testX,
    #                                {"category_output": testCategoryY, "color_output": testColorY}),
    #               epochs=EPOCHS,
    #               verbose=1)
    # # save the model to disk
    # print("[INFO] serializing network...")
    # model.save(args["model"], save_format="h5")
    #
    # # save the category binarizer to disk
    # print("[INFO] serializing category label binarizer...")
    # f = open(args["categorybin"], "wb")
    # f.write(pickle.dumps(categoryLB))
    # f.close()
    # # save the color binarizer to disk
    # print("[INFO] serializing color label binarizer...")
    # f = open(args["colorbin"], "wb")
    # f.write(pickle.dumps(colorLB))
    # f.close()
    #
    # # ---------------------------------------------------------------------------------------
    # # plot the total loss, category loss, and color loss
    # lossNames = ["loss", "category_output_loss", "color_output_loss"]
    # plt.style.use("ggplot")
    # (fig, ax) = plt.subplots(3, 1, figsize=(13, 13))
    #
    # # loop over the loss names
    # for (i, l) in enumerate(lossNames):
    #     # plot the loss for both the training and validation data
    #     title = "Loss for {}".format(l) if l != "loss" else "Total loss"
    #     ax[i].set_title(title)
    #     ax[i].set_xlabel("Epoch #")
    #     ax[i].set_ylabel("Loss")
    #     ax[i].plot(np.arange(0, EPOCHS), H.history[l], label=l)
    #     ax[i].plot(np.arange(0, EPOCHS), H.history["val_" + l],
    #                label="val_" + l)
    #     ax[i].legend()
    #
    # # save the losses figure
    # plt.tight_layout()
    # plt.savefig("{}_losses.png".format(args["plot"]))
    # plt.close()
    #
    # # create a new figure for the accuracies
    # accuracyNames = ["category_output_accuracy", "color_output_accuracy"]
    # plt.style.use("ggplot")
    # (fig, ax) = plt.subplots(2, 1, figsize=(8, 8))
    # # loop over the accuracy names
    # for (i, l) in enumerate(accuracyNames):
    #     # plot the loss for both the training and validation data
    #     ax[i].set_title("Accuracy for {}".format(l))
    #     ax[i].set_xlabel("Epoch #")
    #     ax[i].set_ylabel("Accuracy")
    #     ax[i].plot(np.arange(0, EPOCHS), H.history[l], label=l)
    #     ax[i].plot(np.arange(0, EPOCHS), H.history["val_" + l],
    #                label="val_" + l)
    #     ax[i].legend()
    # # save the accuracies figure
    # plt.tight_layout()
    # plt.savefig("{}_accs.png".format(args["plot"]))
    # plt.close()


if __name__ == "__main__":
    main()









