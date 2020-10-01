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
    Attention

from tensorflow.keras.activations import relu
import keras.backend as K


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


















