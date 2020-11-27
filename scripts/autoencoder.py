#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 11:46:49 2020

@author: labuser
"""



from tensorflow.keras.layers import Conv1D, \
     Softmax, Add, Dense, Activation,\
     Flatten, Reshape

from tensorflow.keras.models import Model
from tensorflow.keras import Input

import keras.backend as K
import numpy as np

   
def autoencoderConv1D(self, signal_shape=(50, 1)):
    """
    Conv2D auto-encoder model.
    Arguments:
        img_shape: e.g. (28, 28, 1) for MNIST
    return:
        (autoencoder, encoder), Model of autoencoder and model of encoder
    """
    input_img = Input(shape=signal_shape)
    # Encoder
    x = Conv1D(48, 3, activation='relu', padding='same')(input_img)
    x = Conv1D(16, 3, activation='relu', padding='same')(x)
    x = Conv1D(16, 3, activation='relu', padding='same')(x)
    shape_before_flattening = K.int_shape(x)
    # at this point the representation is (4, 4, 8) i.e. 128-dimensional
    x = Flatten()(x)
    encoded = Dense(10, activation='relu', name='encoded')(x)

    # Decoder
    x = Dense(np.prod(shape_before_flattening[1:]),
                activation='relu')(encoded)
    # Reshape into an image of the same shape as before our last `Flatten` layer
    x = Reshape(shape_before_flattening[1:])(x)

    x = Conv1D(16, 3, activation='relu', padding='same')(x)
    #x = UpSampling1D((2, 2))(x)
    x = Conv1D(16, 3, activation='relu', padding='same')(x)
    #x = UpSampling1D((2, 2))(x)
    x = Conv1D(48, 3, activation='relu', padding='same')(x)
    #x = UpSampling1D((2, 2))(x)
    decoded = Conv1D(1, 3, activation='sigmoid', padding='same')(x)

    return Model(inputs=input_img, outputs=decoded, name='AE'), Model(inputs=input_img, outputs=encoded, name='encoder')
