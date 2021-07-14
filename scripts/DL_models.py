#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 10:49:04 2020

@author: labuser
"""

"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Deepbinner/

This file is part of Deepbinner. Deepbinner is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Deepbinner is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Deepbinner.
If not, see <http://www.gnu.org/licenses/>.
"""

from tensorflow.keras.layers import Conv1D, MaxPooling1D, AveragePooling1D, Dropout, concatenate,\
    BatchNormalization, GaussianNoise, GlobalAveragePooling1D, Softmax, Add, Dense, Activation,\
    Attention, Flatten, LayerNormalization

from tensorflow.keras.activations import relu
import keras.backend as K
import tensorflow as tf


def build_deepbinner(inputs, class_count):
    """
    This function builds the standard network used in Deepbinner.
    """
    x = inputs

    # Add some noise to augment the training data.
    #x = GaussianNoise(stddev=0.02)(x)

    # Conv layer with stride of 2 (halves the size)
    x = Conv1D(filters=48, kernel_size=3, strides=2, padding='same', activation='relu')(x)

    x = BatchNormalization()(x)
    x = Dropout(rate=0.15)(x)

    # Conv group: 3 layers of 3-kernels
    x = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x)
    x = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x)
    x = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x)
    x = MaxPooling1D(pool_size=2)(x)

    x = BatchNormalization()(x)
    x = Dropout(rate=0.15)(x)

    # Bottleneck down to 16 filters (reduces the number of parameters a bit)
    x = Conv1D(filters=16, kernel_size=1, activation='relu')(x)

    # Conv group: 2 layers of 3-kernels
    x = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x)
    x = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x)
    x = MaxPooling1D(pool_size=2)(x)

    x = BatchNormalization()(x)
    x = Dropout(rate=0.15)(x)

    # Conv group: 2 layers of 3-kernels
    x = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x)
    x = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x)
    x = MaxPooling1D(pool_size=2)(x)

    x = BatchNormalization()(x)
    x = Dropout(rate=0.15)(x)

    # Inception-style group
    x1 = AveragePooling1D(pool_size=3, strides=1, padding='same')(x)
    x1 = Conv1D(filters=48, kernel_size=1, padding='same', activation='relu')(x1)
    x2 = Conv1D(filters=48, kernel_size=1, padding='same', activation='relu')(x)
    x3 = Conv1D(filters=16, kernel_size=1, padding='same', activation='relu')(x)
    x3 = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x3)
    x4 = Conv1D(filters=16, kernel_size=1, padding='same', activation='relu')(x)
    x4 = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x4)
    x4 = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x4)
    x = concatenate([x1, x2, x3, x4], axis=2)
    x = MaxPooling1D(pool_size=2)(x)

    x = BatchNormalization()(x)
    x = Dropout(rate=0.15)(x)

    # Conv layer with stride of 2 (halves the size)
    x = Conv1D(filters=48, kernel_size=3, strides=2, activation='relu', padding='same')(x)

    x = BatchNormalization()(x)
    x = Dropout(rate=0.15)(x)

    # Conv group: 2 layers of 3-kernels
    x = Conv1D(filters=48, kernel_size=3, activation='relu', padding='same')(x)
    x = Conv1D(filters=48, kernel_size=3, activation='relu', padding='same')(x)
    x = MaxPooling1D(pool_size=2)(x)

    x = BatchNormalization()(x)
    x = Dropout(rate=0.15)(x)

    # Finish with a global average pooling approach (no fully connected layers)
    x = Conv1D(filters=class_count, kernel_size=1, activation='relu')(x)
    x = GlobalAveragePooling1D()(x)
    x = Dense(1, activation='sigmoid')(x)

    return x

def _bn_relu(input):
    """Helper to build a BN -> relu block
    """
    norm = BatchNormalization()(input)
    return Activation("relu")(norm)

def build_Jasper(inputs, Deep=None):
    '''
    a Jasper BxR model has B blocks, each with R subblocks. 
    Each sub-block applies the following operations: 
        1Dconvolution, 
        batch norm, 
        ReLU, 
        dropout.
    All sub-blocks in a block have the same number of output channels.
    Each block input is connected directly into the last subblock via 
    a residual connection. The residual connection is first projected 
    through a 1x1 convolution to account for different numbers of input 
    and output channels, then through a batch norm layer. The output
    of this batch norm layer is added to the output of the batch norm
    layer in the last sub-block. The result of this sum is passed 
    through the activation function and dropout to produce the output 
    of the current block.
    All Jasper models have four additional convolutional
    blocks: one pre-processing and three post-processing.
    https://arxiv.org/pdf/1904.03288.pdf
    '''
    x = inputs
    '''
    # pre-processing (prolog) Conv layer
    x = Conv1D(filters=48, 
               kernel_size=3,
               padding='same',
               dilation_rate=2
               )(x)
    x = _bn_relu(x)
    '''
    x1 = Conv1D(filters=48, kernel_size=1, padding='same', activation='relu')(x)
    x2 = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x)
    x3 = Conv1D(filters=48, kernel_size=3, dilation_rate=2, padding='same', activation='relu')(x)
    x = concatenate([x1, x2, x3], axis=2)
    x = MaxPooling1D(pool_size=3)(x)
    x = _bn_relu(x)
    
    ## First block
    
    n_filter = 48
    kernel_s = 2
    droppout = 0.15
    
    a = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
              )(x)
    a = _bn_relu(a)
    a = Dropout(rate=droppout)(a)
    
    a = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(a)
    a = _bn_relu(a)
    a = Dropout(rate=droppout)(a)
    
    a = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(a)
    a = BatchNormalization()(a)
    '''
    #  1X1 conv residual connection
    '''
    x1 = Conv1D(filters=48, 
               kernel_size=2, 
               padding='same', 
               )(x)
    x1 = BatchNormalization()(x1)
    '''
    #  1X1 conv residual connection
    '''
    a = Add()([x1, a])
    a = Activation("relu")(a)
    a = Dropout(rate=droppout)(a)
    #a = MaxPooling1D(pool_size=2,padding='same')(a)

    ## Second block
    n_filter = 48
    kernel_s = 3
    droppout = 0.15
    
    b = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(a)
    b = BatchNormalization()(b)
    b = Activation("relu")(b)
    b = Dropout(rate=droppout)(b)
    
    b = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(b)
    b =_bn_relu(b)
    b = Dropout(rate=droppout)(b)
    
    b = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(b)
    b = BatchNormalization()(b)

    '''
    #  1X1 conv residual connection
    '''
    x1 = Conv1D(filters=48, 
               kernel_size=1, 
               padding='same', 
               )(a)
    x1 = BatchNormalization()(x1)
    '''
    #  1X1 conv residual connection
    '''
    if Deep == True:
        b = Add()([x1, b, a])
    else:
        b = Add()([x1, b])
    b = Activation("relu")(b)
    b = Dropout(rate=droppout)(b)
    #b = MaxPooling1D(pool_size=2,padding='same')(b)

    ## Third block
    n_filter = 48
    kernel_s = 3
    droppout = 0.15
    
    c = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(b)
    c = _bn_relu(c)
    c = Dropout(rate=droppout)(c)
    
    c = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(c)
    c = _bn_relu(c)
    c = Dropout(rate=droppout)(c)

    c = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(c)
    c = BatchNormalization()(c)

    '''
    #  1X1 conv residual connection
    '''
    x1 = Conv1D(filters=48, 
               kernel_size=1, 
               padding='same', 
               )(b)
    x1 = BatchNormalization()(x1)
    '''
    #  1X1 conv residual connection
    '''
    if Deep == True:
        c = Add()([x1, b, a, c])
    else:
        c = Add()([x1, c])
    c = Activation("relu")(c)
    c = Dropout(rate=droppout)(c)
    #c = MaxPooling1D(pool_size=2,padding='same')(c)

    ## fourth block
    n_filter = 48
    kernel_s = 3
    droppout = 0.15
    
    d = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(c)
    d =_bn_relu(d)
    d = Dropout(rate=droppout)(d) 
    
    d = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(d)
    d =_bn_relu(d)
    d = Dropout(rate=droppout)(d)

    d = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(d)
    d = BatchNormalization()(d)

    '''
    #  1X1 conv residual connection
    '''
    x1 = Conv1D(filters=48, 
               kernel_size=2, 
               padding='same', 
               )(c)
    x1 = BatchNormalization()(x1)
    '''
    #  1X1 conv residual connection
    '''
    if Deep == True:
        d = Add()([x1, b, a, c, d])
    else:
        d = Add()([x1, d])
    d = Activation("relu")(d)
    d = Dropout(rate=droppout)(d)
    #d = MaxPooling1D(pool_size=2,padding='same')(d)

    # epilog conv layer
    x = Conv1D(filters=48, 
               kernel_size=3,
               dilation_rate=2,
               padding='same')(d)
    x =_bn_relu(x)
    
    x = Conv1D(filters=48, 
               kernel_size=3,
               dilation_rate=2,
               padding='same')(x)
    x =_bn_relu(x)
    
    x = Conv1D(filters=1, 
               kernel_size=1, 
               padding='same', 
               activation='relu')(x)
    
    x = GlobalAveragePooling1D()(x)
    x = Dense(1, activation='sigmoid')(x)
    
    return x

def build_Jasper_2inputs(input1, input2, Deep=None):
    '''
    a Jasper BxR model has B blocks, each with R subblocks. 
    Each sub-block applies the following operations: 
        1Dconvolution, 
        batch norm, 
        ReLU, 
        dropout.
    All sub-blocks in a block have the same number of output channels.
    Each block input is connected directly into the last subblock via 
    a residual connection. The residual connection is first projected 
    through a 1x1 convolution to account for different numbers of input 
    and output channels, then through a batch norm layer. The output
    of this batch norm layer is added to the output of the batch norm
    layer in the last sub-block. The result of this sum is passed 
    through the activation function and dropout to produce the output 
    of the current block.
    All Jasper models have four additional convolutional
    blocks: one pre-processing and three post-processing.
    https://arxiv.org/pdf/1904.03288.pdf
    '''
    x = input1
    '''
    # pre-processing (prolog) Conv layer
    x = Conv1D(filters=48, 
               kernel_size=3,
               padding='same',
               dilation_rate=2
               )(x)
    x = _bn_relu(x)
    '''
    x1 = Conv1D(filters=48, kernel_size=1, padding='same', activation='relu')(x)
    x2 = Conv1D(filters=48, kernel_size=3, padding='same', activation='relu')(x)
    x3 = Conv1D(filters=48, kernel_size=3, dilation_rate=2, padding='same', activation='relu')(x)
    x = concatenate([x1, x2, x3], axis=2)
    x = MaxPooling1D(pool_size=3)(x)
    x = _bn_relu(x)
    
    I1 = Conv1D(filters=12, kernel_size=1, padding='same', activation='relu')(input2)
    I2 = Conv1D(filters=12, kernel_size=2, padding='same', activation='relu')(input2)
    I3 = Conv1D(filters=12, kernel_size=2, padding='same', activation='relu')(input2)
    I = concatenate([I1, I2, I3], axis=2)
    I = _bn_relu(I)
    
    ## First block
    
    n_filter = 48
    kernel_s = 2
    droppout = 0.15
    
    a = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
              )(x)
    a = _bn_relu(a)
    a = Dropout(rate=droppout)(a)
    
    a = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(a)
    a = _bn_relu(a)
    a = Dropout(rate=droppout)(a)
    
    a = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(a)
    a = BatchNormalization()(a)
    '''
    #  1X1 conv residual connection
    '''
    x1 = Conv1D(filters=48, 
               kernel_size=2, 
               padding='same', 
               )(x)
    x1 = BatchNormalization()(x1)
    '''
    #  1X1 conv residual connection
    '''
    a = Add()([x1, a])
    a = Activation("relu")(a)
    a = Dropout(rate=droppout)(a)
    #a = MaxPooling1D(pool_size=2,padding='same')(a)

    ## Second block
    n_filter = 48
    kernel_s = 3
    droppout = 0.15
    
    b = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(a)
    b = BatchNormalization()(b)
    b = Activation("relu")(b)
    b = Dropout(rate=droppout)(b)
    
    b = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(b)
    b =_bn_relu(b)
    b = Dropout(rate=droppout)(b)
    
    b = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(b)
    b = BatchNormalization()(b)

    '''
    #  1X1 conv residual connection
    '''
    x1 = Conv1D(filters=48, 
               kernel_size=1, 
               padding='same', 
               )(a)
    x1 = BatchNormalization()(x1)
    '''
    #  1X1 conv residual connection
    '''
    if Deep == True:
        b = Add()([x1, b, a])
    else:
        b = Add()([x1, b])
    b = Activation("relu")(b)
    b = Dropout(rate=droppout)(b)
    #b = MaxPooling1D(pool_size=2,padding='same')(b)

    ## Third block
    n_filter = 48
    kernel_s = 3
    droppout = 0.15
    
    c = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(b)
    c = _bn_relu(c)
    c = Dropout(rate=droppout)(c)
    
    c = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(c)
    c = _bn_relu(c)
    c = Dropout(rate=droppout)(c)

    c = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(c)
    c = BatchNormalization()(c)

    '''
    #  1X1 conv residual connection
    '''
    x1 = Conv1D(filters=48, 
               kernel_size=1, 
               padding='same', 
               )(b)
    x1 = BatchNormalization()(x1)
    '''
    #  1X1 conv residual connection
    '''
    if Deep == True:
        c = Add()([x1, b, a, c])
    else:
        c = Add()([x1, c])
    c = Activation("relu")(c)
    c = Dropout(rate=droppout)(c)
    #c = MaxPooling1D(pool_size=2,padding='same')(c)

    ## fourth block
    n_filter = 48
    kernel_s = 3
    droppout = 0.15
    
    d = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(c)
    d =_bn_relu(d)
    d = Dropout(rate=droppout)(d) 
    
    d = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(d)
    d =_bn_relu(d)
    d = Dropout(rate=droppout)(d)

    d = Conv1D(filters=n_filter, 
               kernel_size=kernel_s,
               padding='same',
               )(d)
    d = BatchNormalization()(d)

    '''
    #  1X1 conv residual connection
    '''
    x1 = Conv1D(filters=48, 
               kernel_size=2, 
               padding='same', 
               )(c)
    x1 = BatchNormalization()(x1)
    '''
    #  1X1 conv residual connection
    '''
    if Deep == True:
        d = Add()([x1, b, a, c, d])
    else:
        d = Add()([x1, d])
    d = Activation("relu")(d)
    d = Dropout(rate=droppout)(d)
    #d = MaxPooling1D(pool_size=2,padding='same')(d)

    # epilog conv layer
    x = Conv1D(filters=48, 
               kernel_size=3,
               dilation_rate=2,
               padding='same')(d)
    x =_bn_relu(x)
    
    x = Conv1D(filters=48, 
               kernel_size=3,
               dilation_rate=2,
               padding='same')(x)
    x =_bn_relu(x)
    
    x = Conv1D(filters=1, 
               kernel_size=1, 
               padding='same', 
               activation='relu')(x)
    
    x = concatenate([Flatten()(x), Flatten()(I)])
    
    x = Dense(25)(x)
    
    x = Dense(1, activation='sigmoid')(x)
    
    return x



    
    
def build_TCN(inputs, output):
    '''
    '''
    x1 = TCN(return_sequences=False,
             nb_filters=48,
             nb_stacks=4,
             kernel_size=3,
             dilations=[1,2,4,8], 
             padding='same',
             activation='relu',
             use_batch_norm=True,
             dropout_rate=0.2)(inputs)
    
    output = Dense(output, activation='sigmoid')(x1)
    
    return output
    

def build_TCN_cc(inputs, output):
    '''
    '''
    x1 = TCN(nb_filters=48,
             kernel_size=2,
             dilations=[1], 
             padding='same',
             use_batch_norm=True,
             return_sequences=True,
             dropout_rate=0.2)(inputs)
     
    x1 = TCN(return_sequences=True,
             nb_filters=48,
             kernel_size=2,
             dilations=[1], 
             padding='same',
             use_batch_norm=True,
             dropout_rate=0.2)(x1)
    
    x1 = TCN(return_sequences=True,
             nb_filters=48,
             kernel_size=2,
             dilations=[1], 
             padding='same',
             use_batch_norm=True,
             dropout_rate=0.2)(x1)
    '''
    residual connection
    '''
    xres = Conv1D(filters=48, 
                kernel_size=1, 
                padding='same', 
                )(inputs)
    '''
    residual connection
    '''
    
    x1 = Add()([x1, xres])
    x1 = Activation("relu")(x1)
    
    x2 = TCN(nb_filters=48,
             kernel_size=2,
             dilations=[2], 
             padding='same',
             use_batch_norm=True,
             return_sequences=True,
             dropout_rate=0.2)(x1)
     
    x2 = TCN(return_sequences=True,
             nb_filters=48,
             kernel_size=2,
             dilations=[2], 
             padding='same',
             use_batch_norm=True,
             dropout_rate=0.2)(x2)
    
    x2 = TCN(return_sequences=True,
             nb_filters=48,
             kernel_size=2,
             dilations=[2], 
             padding='same',
             use_batch_norm=True,
             dropout_rate=0.2)(x2)
    '''
    residual connection
    '''
    xres = Conv1D(filters=48, 
                kernel_size=1, 
                padding='same', 
                )(x1)
    '''
    residual connection
    '''
    
    x2 = Add()([x2, xres])
    x2 = Activation("relu")(x2)
    
    x3 = TCN(nb_filters=48,
             kernel_size=2,
             dilations=[4], 
             padding='same',
             use_batch_norm=True,
             return_sequences=True,
             dropout_rate=0.2)(x2)
     
    x3 = TCN(return_sequences=True,
             nb_filters=48,
             kernel_size=2,
             dilations=[4], 
             padding='same',
             use_batch_norm=True,
             dropout_rate=0.2)(x3)
    
    x3 = TCN(return_sequences=True,
             nb_filters=48,
             kernel_size=2,
             dilations=[4], 
             padding='same',
             use_batch_norm=True,
             dropout_rate=0.2)(x3)
    '''
    residual connection
    '''
    xres = Conv1D(filters=48, 
                kernel_size=1, 
                padding='same', 
                )(x2)
    '''
    residual connection
    '''
    
    x3 = Add()([x3, xres])
    x3 = Activation("relu")(x3)
    
    x4 = TCN(nb_filters=48,
             kernel_size=2,
             dilations=[8], 
             padding='same',
             use_batch_norm=True,
             return_sequences=True,
             dropout_rate=0.2)(x3)
     
    x4 = TCN(return_sequences=True,
             nb_filters=48,
             kernel_size=2,
             dilations=[8], 
             padding='same',
             use_batch_norm=True,
             dropout_rate=0.2)(x4)
    
    x4 = TCN(return_sequences=True,
             nb_filters=48,
             kernel_size=2,
             dilations=[8], 
             padding='same',
             use_batch_norm=True,
             dropout_rate=0.2)(x4)
    '''
    residual connection
    '''
    xres = Conv1D(filters=48, 
                kernel_size=1, 
                padding='same', 
                )(x3)
    '''
    residual connection
    '''
    x4 = Add()([x4, xres])
    x4= Activation("relu")(x4)
    
    x4 = GlobalAveragePooling1D()(x4)
    
    output = Dense(output, activation='sigmoid')(x4)
    
    return output


from tensorflow.keras.layers import  Reshape

from tensorflow.keras.models import Model
from tensorflow.keras import Input

import numpy as np

def autoencoderConv1D(signal_shape=(50, 1)):
    """
    Conv2D auto-encoder model.
    Arguments:
        img_shape: e.g. (28, 28, 1) for MNIST
    return:
        (autoencoder, encoder), Model of autoencoder and model of encoder
    """
    input_img = Input(shape=signal_shape)
    # Encoder
    x = Conv1D(32, 3, activation='relu', padding='same')(input_img)
    x = Conv1D(24, 3, activation='relu', padding='same')(x)
    x = Conv1D(16, 3, activation='relu', padding='same')(x)
    shape_before_flattening = K.int_shape(x)
    # at this point the representation is (4, 4, 8) i.e. 128-dimensional
    x = Flatten()(x)
    encoded = Dense(50, activation='relu', name='encoded')(x)

    # Decoder
    x = Dense(np.prod(shape_before_flattening[1:]),
                activation='relu')(encoded)
    # Reshape into an image of the same shape as before our last `Flatten` layer
    x = Reshape(shape_before_flattening[1:])(x)
    
    x = Conv1D(16, 3, activation='relu', padding='same')(x)
    #x = UpSampling1D((2, 2))(x)
    x = Conv1D(24, 3, activation='relu', padding='same')(x)
    #x = UpSampling1D((2, 2))(x)
    x = Conv1D(32, 3, activation='relu', padding='same')(x)
    #x = UpSampling1D((2, 2))(x)
    decoded = Conv1D(1, 2, activation='relu', padding='same')(x)

    return Model(inputs=input_img, outputs=decoded, name='AE'), \
           Model(inputs=input_img, outputs=encoded, name='encoder')














