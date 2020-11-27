#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 09:08:56 2020

@author: labuser
"""

import pandas as pd
import numpy as np

# Make numpy values easier to read.
np.set_printoptions(precision=3, suppress=True)

import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras.layers.experimental import preprocessing
import os

train_dataset_url = "https://storage.googleapis.com/download.tensorflow.org/data/iris_training.csv"

train_dataset_fp = tf.keras.utils.get_file(fname=os.path.basename(train_dataset_url),

                                           origin=train_dataset_url)
batch_size = 32
column_names = ['sepal_length', 'sepal_width', 'petal_length', 'petal_width', 'species']

feature_names = column_names[:-1]
label_name = column_names[-1]


train_dataset = tf.data.experimental.make_csv_dataset(
                                train_dataset_fp,
                                batch_size,
                                column_names=column_names,
                                label_name=label_name,
                                num_epochs=1)



def pack_features_vector(features, labels):
  """Pack the features into a single array."""
  features = tf.stack(list(features.values()), axis=1)
  return features, labels

path_signal = '/media/labuser/Data/nanopore/m6A_classifier/data/Epinano/doubleAA/TRAIN.csv'

train_dataset = tf.data.experimental.make_csv_dataset(
                                                    path_signal,
                                                    batch_size=1000,
                                                    label_name='label',
                                                    num_epochs=1,
                                                    select_columns=['event', 'distances', 'dwell', 'label'],
                                                    column_defaults= ['float','float','float','int64'],
                                                    field_delim='\t')

train_dataset = train_dataset.map(pack_features_vector)

#df = pd.read_csv(path_signal, sep='\t')


for x, y in train_dataset:
    print(x)
    print(np.array(x).reshape(10, 100, 3))
    print(y)
    break

 # Define model Inout
inputs = Input(shape=(None, 3))

output = build_Jasper(inputs,Deep=True)
#output = build_deepbinner(inputs, 1)

model = Model(inputs=inputs, outputs=output)

model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=['accuracy',
                       keras_metrics.recall(),
                       keras_metrics.precision(),
                       ])

#checkpoint = ModelCheckpoint('/media/labuser/Data/nanopore/m6A_classifier/scripts/models/m6A_classifier_DB.h5',
#                             monitor='val_acc', verbose=1,
#                             save_best_only=True, mode='max')

OUT_FOLDER = '/media/labuser/Data/nanopore/m6A_classifier/plots/'

model.summary()

epochs =1
histories = []

for e in epoch:
    for x, y in train_dataset:
        #Train the model on batches
        
        history = model.fit(x, y)
        
    
        
        
zip_path = tf.keras.utils.get_file(
origin='https://storage.googleapis.com/download.tensorflow.org/data/creditcard.zip',
fname='creditcard.zip',
extract=True)
csv_path = zip_path.replace('.zip', '.csv')

creditcard_ds = tf.data.experimental.make_csv_dataset(
csv_path, batch_size=1024, label_name="Class",
# Set the column types: 30 floats and an int.
column_defaults=[float()]*30+[int()])

csv_path = zip_path.replace('.zip', '.csv')

        
        
        
            
            
            
            
            
            
            
            
            
            

