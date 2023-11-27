# created: 03-23-23
# author: PD 
# imports
import glob
import os
import sys
sys.path.append('/content/drive/MyDrive/Spindles/')

from scipy.io import savemat,loadmat
import numpy as np

import tensorflow as tf
from tensorflow.python.ops.math_ops import xlogy
from tensorflow.keras.layers import Input, Dense, Dropout
from tensorflow.keras.models import Model
from keras.constraints import UnitNorm, Constraint
from keras.models import Model

#from processp import matchSpiRip


## MODEL CREATION


# Common method to tie encoder and decoder layers, so that they are symmetric
class TiedDense(tf.keras.layers.Layer):
    def __init__(self, tied_layer, activation=None, **kwargs):
        self.tied_layer = tied_layer
        self.activation = tf.keras.activations.get(activation)
        super(TiedDense, self).__init__(**kwargs)

    def build(self, input_shape):
        self.bias = self.add_weight(name="bias", shape=self.tied_layer.input_shape[-1], initializer=self.tied_layer.bias_initializer, trainable=True)
        #self.set_weights(self.transpose_weights)
        super(TiedDense, self).build(input_shape)

    def call(self, inputs):
        output = tf.matmul(inputs, self.tied_layer.weights[0], transpose_b=True)
        if self.activation is not None:
            output = self.activation(output+self.bias)
        return output

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.tied_layer.units)




def SimpleAutoencoder( code_size, inp, acti=None,ls=[]):


    x=inp

    #ls=[126,140,100,code_size]
    enc=[]
    for k,i in enumerate(ls[1:]):
        enc.append(Dense(i,activation=acti,name="enc{}".format(k),kernel_constraint=UnitNorm(axis=0)))

    dec=[]
    for k,i in enumerate(enc[::-1]):
        dec.append(TiedDense(i,activation=acti,name="dec{}".format(k)))

    for enci in enc:
        x=enci(x)
    encode=x
    for deci in dec:
        x=deci(x)


    #extra conditioning to seperate global and local spindles
    d2 = Dense(2,activation="softmax",name="class")(encode)

    model = Model(inputs=inp, outputs=[x,d2])
    model.summary()
    return model

recon_loss = 'mse'
class_loss = 'categorical_crossentropy'

metrics = ['accuracy']

input_shape = (126,)
code_size = 10
inp = Input(shape=input_shape)

#initialising the model
ls=[126,100,code_size]
model = SimpleAutoencoder(code_size, inp=inp,ls=ls,acti="tanh")

#using mse loss for autoencoder and binary_crossentropy for global/local classification task
model.compile(optimizer=tf.keras.optimizers.Adam(1e-3), loss=["mse","binary_crossentropy"],loss_weights=[1000,0.1],metrics=[metrics,metrics])


## DATA PROCESSING


# Create a list of all data file paths
data_path = "/content/drive/MyDrive/Spindles/Raw_Data/forPID/OR17/Hilbert_norm/*.npz"
data_files = glob.glob(data_path)

# Define a function to load data from a single file
def load_data(file_path):
    with np.load(file_path) as data:
        x = data["x"]
        y = data["y"]
    return x, y

x_data = np.empty((0, 126), dtype=np.float64)
y_data = np.empty((0, 1), dtype=np.float64)


# Load data from each file and create a dataset of (x, y) pairs
for file_path in data_files:
    x, y = load_data(file_path)
    x_data = np.concatenate([x_data, x], axis=0)
    y_data = np.concatenate([y_data, y[:,None]], axis=0)


from sklearn.preprocessing import OneHotEncoder
ohe=OneHotEncoder(sparse=False)

#converting 
y_dataq =ohe.fit_transform(y_data)

train_log_path = "model2022_train_log.csv"
model_path = "model2022.h5"

#train model with batch size=64 for 100epochs
model.fit(x_data,[x_data,y_dataq],epochs=100,batch_size=64,validation_split=0.1,

        callbacks=[
                tf.keras.callbacks.ModelCheckpoint(model_path, save_best_only=True,save_weights_only=True),
                tf.keras.callbacks.CSVLogger(train_log_path)],

          )


## INFERENCE AND ANALYSIS

model.load_weights(model_path)

#Load only the encoder 
encoder = Model(inputs=inp, outputs=model.layers[-5].output)
#compress the data 
encoded=encoder(x_data)

#Fitting PCA for visualisation
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

pca = PCA(n_components=10)
pca_data = pca.fit_transform(encoded)

pca1 = pca_data[:, 0]
pca2 = pca_data[:, 1]
plt.scatter(pca1, pca2,c=y_data)
plt.xlabel('PCA 1')
plt.ylabel('PCA 2')
plt.show()
