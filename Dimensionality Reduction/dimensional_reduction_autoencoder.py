from __future__ import division, print_function, absolute_import

import time

from util import read_csv, save_csv_data, draw_loss_function
import random

import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import csv

# Import MNIST data
# from tensorflow.examples.tutorials.mnist import input_data

# mnist = input_data.read_data_sets("/tmp/data/", one_hot=True)

target = 'flowcap'

# IMPC_data = '/Local-Scratch/CMPT-884-workspace/flowCAP-II/result/feat/file-cell-countAdj.csv'
if target == 'IMPC':
    IMPC_dir = '/Local-Scratch/CMPT-884-workspace/IMPC/result/P1/Sanger_SPLEEN/feat/file-cell-countAdj.csv'
    data, sample_number_list = read_csv(IMPC_dir)
elif target == 'flowcap':
    flow_cap_dir = '/Local-Scratch/CMPT-884-workspace/flowCAP-II/result/feat/file-cell-countAdj.csv'
    data, sample_number_list = read_csv(flow_cap_dir)
else:
    raise ValueError('wrong target')

shape = data.shape
feature_num = shape[1]
data_num = shape[0]

# Training Parameters
learning_rate_list = [0.0001, 0.001, 0.01]
iteration_number = 100
batch_size = 100

display_step = 1000
examples_to_show = 10

# Network Parameters
num_hidden_1 = int(feature_num / 4)  # 1st layer num features
num_hidden_2 = int(feature_num / 8)  # 2nd layer num features (the latent dim)
num_input = feature_num  # MNIST data input (img shape: 28*28)

# tf Graph input (only pictures)
X = tf.placeholder("float", [None, num_input])

weights = {
    'encoder_h1': tf.Variable(tf.random_normal([num_input, num_hidden_1])),
    'encoder_h2': tf.Variable(tf.random_normal([num_hidden_1, num_hidden_2])),
    'decoder_h1': tf.Variable(tf.random_normal([num_hidden_2, num_hidden_1])),
    'decoder_h2': tf.Variable(tf.random_normal([num_hidden_1, num_input])),
}
biases = {
    'encoder_b1': tf.Variable(tf.random_normal([num_hidden_1])),
    'encoder_b2': tf.Variable(tf.random_normal([num_hidden_2])),
    'decoder_b1': tf.Variable(tf.random_normal([num_hidden_1])),
    'decoder_b2': tf.Variable(tf.random_normal([num_input])),
}


# Building the encoder
def encoder(x):
    # Encoder Hidden layer with sigmoid activation #1
    layer_1 = tf.nn.sigmoid(tf.add(tf.matmul(x, weights['encoder_h1']),
                                   biases['encoder_b1']))
    # Encoder Hidden layer with sigmoid activation #2
    layer_2 = tf.nn.sigmoid(tf.add(tf.matmul(layer_1, weights['encoder_h2']),
                                   biases['encoder_b2']))
    return layer_2


# Building the decoder
def decoder(x):
    # Decoder Hidden layer with sigmoid activation #1
    layer_1 = tf.nn.sigmoid(tf.add(tf.matmul(x, weights['decoder_h1']),
                                   biases['decoder_b1']))
    # Decoder Hidden layer with sigmoid activation #2
    layer_2 = tf.nn.sigmoid(tf.add(tf.matmul(layer_1, weights['decoder_h2']),
                                   biases['decoder_b2']))
    return layer_2


plt.figure(figsize=(7, 5))
for learning_rate in learning_rate_list:
    start_time = time.time()
    # Construct model
    encoder_op = encoder(X)
    decoder_op = decoder(encoder_op)

    # Prediction
    y_pred = decoder_op
    # Targets (Labels) are the input data.
    y_true = X

    # Define loss and optimizer, minimize the squared error
    loss = tf.sqrt(tf.reduce_mean(tf.pow(y_true - y_pred, 2)))
    optimizer = tf.train.AdamOptimizer(learning_rate).minimize(loss)

    # Initialize the variables (i.e. assign their default value)
    init = tf.global_variables_initializer()

    # Start Training
    # Start a new TF session
    with tf.Session() as sess:
        # Run the initializer
        sess.run(init)

        # Training

        iteration_number_plot = []
        loss_plot = []

        for i in range(1, iteration_number + 1):
            training_indexs = range(0, data_num)
            random.shuffle(training_indexs)
            splits_number = int(data_num / batch_size)

            for number in range(0, splits_number + 1):
                if number == splits_number:
                    batch_indexs = training_indexs[number * batch_size:]
                else:
                    batch_indexs = training_indexs[number * batch_size: (number + 1) * batch_size]

                # Prepare Data
                batch_data = data[batch_indexs]
                # batch_x, _ = mnist.train.next_batch(batch_size)

                # Run optimization op (backprop) and cost op (to get loss value)
                _, l = sess.run([optimizer, loss], feed_dict={X: batch_data})
                # Display logs per step
                # if i % display_step == 0 or i == 1:

            _, l = sess.run([optimizer, loss], feed_dict={X: data})
            loss_plot.append(l)
            iteration_number_plot.append(i)
            print('iteration %i: Minibatch Loss: %f' % (i, l))

        dimension_reduction = (sess.run([encoder_op], feed_dict={X: data}))[0]

    plt.plot(iteration_number_plot, loss_plot, linewidth=3, label='LR={0}'.format(learning_rate))
    print("---sdae for {0} {1} seconds ---".format(target, time.time() - start_time))
plt.xlabel('Iteration', fontsize=18)
plt.ylabel('Loss', fontsize=18)
plt.legend(loc="best",fontsize=15)
plt.savefig('./draw_result/{1}-ae-loss-{0}-lr.png'.format(iteration_number, target))

# save_csv_data(dir='./dr_results/{1}-ae-{0}.csv'.format(int(feature_num / 8), target), data=dimension_reduction,
#               sample_number_list=sample_number_list,
#               name='ae')

# draw_loss_function(iteration_number_plot, loss_plot,
#                    './draw_result/{1}-ae-loss-{0}'.format(iteration_number, target))

print('still working')

# # Testing
# # Encode and decode images from test set and visualize their reconstruction.
# n = 4
# canvas_orig = np.empty((28 * n, 28 * n))
# canvas_recon = np.empty((28 * n, 28 * n))
# for i in range(n):
#     # MNIST test set
#     batch_x, _ = mnist.test.next_batch(n)
#     # Encode and decode the digit image
#     g = sess.run(decoder_op, feed_dict={X: batch_x})
#
#     # Display original images
#     for j in range(n):
#         # Draw the original digits
#         canvas_orig[i * 28:(i + 1) * 28, j * 28:(j + 1) * 28] = \
#             batch_x[j].reshape([28, 28])
#     # Display reconstructed images
#     for j in range(n):
#         # Draw the reconstructed digits
#         canvas_recon[i * 28:(i + 1) * 28, j * 28:(j + 1) * 28] = \
#             g[j].reshape([28, 28])
#
# print("Original Images")
# plt.figure(figsize=(n, n))
# plt.imshow(canvas_orig, origin="upper", cmap="gray")
# plt.show()
#
# print("Reconstructed Images")
# plt.figure(figsize=(n, n))
# plt.imshow(canvas_recon, origin="upper", cmap="gray")
# plt.show()
