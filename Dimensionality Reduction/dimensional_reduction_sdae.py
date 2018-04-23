"""Tutorial on how to create a denoising autoencoder w/ Tensorflow.

Parag K. Mital, Jan 2016
"""
import csv
import random
import time

import tensorflow as tf
import numpy as np
import math
from util import read_csv, save_csv_data, draw_loss_function


# from libs.utils import corrupt


def corrupt(x):
    """Take an input tensor and add uniform masking.
    Parameters
    ----------
    x : Tensor/Placeholder
        Input to corrupt.
    Returns
    -------
    x_corrupted : Tensor
        50 pct of values corrupted.
    """
    return tf.multiply(x, tf.cast(tf.random_uniform(shape=tf.shape(x),
                                                    minval=0,
                                                    maxval=2,
                                                    dtype=tf.int32), tf.float32))


# %%
def autoencoder(dimensions):
    """Build a deep denoising autoencoder w/ tied weights.

    Parameters
    ----------
    dimensions : list, optional
        The number of neurons for each layer of the autoencoder.

    Returns
    -------
    x : Tensor
        Input placeholder to the network
    z : Tensor
        Inner-most latent representation
    y : Tensor
        Output reconstruction of the input
    cost : Tensor
        Overall cost to use for training
    """
    # input to the network
    x = tf.placeholder(tf.float32, [None, dimensions[0]], name='x')

    # Probability that we will corrupt input.
    # This is the essence of the denoising autoencoder, and is pretty
    # basic.  We'll feed forward a noisy input, allowing our network
    # to generalize better, possibly, to occlusions of what we're
    # really interested in.  But to measure accuracy, we'll still
    # enforce a training signal which measures the original image's
    # reconstruction cost.
    #
    # We'll change this to 1 during training
    # but when we're ready for testing/production ready environments,
    # we'll put it back to 0.
    corrupt_prob = tf.placeholder(tf.float32, [1])
    current_input = corrupt(x) * corrupt_prob + x * (1 - corrupt_prob)

    # Build the encoder
    # encoder = []
    for layer_i, n_output in enumerate(dimensions[1:]):
        n_input = int(current_input.get_shape()[1])
        W = tf.Variable(
            tf.random_normal([n_input, n_output]))
        b = tf.Variable(tf.zeros([n_output]))
        # encoder.append(W)
        output = tf.nn.tanh(tf.matmul(current_input, W) + b)
        current_input = output
    # latent representation
    z = current_input
    # encoder.reverse()
    # Build the decoder using the same weights
    for layer_i, n_output in enumerate(dimensions[:-1][::-1]):
        # W = tf.transpose(encoder[layer_i])
        n_input = int(current_input.get_shape()[1])
        W = tf.Variable(
            tf.random_normal([n_input, n_output]))
        b = tf.Variable(tf.zeros([n_output]))
        output = tf.nn.tanh(tf.matmul(current_input, W) + b)
        current_input = output
    # now have the reconstruction through the network
    y = current_input
    # cost function measures pixel-wise difference
    cost = tf.sqrt(tf.reduce_mean(tf.square(y - x)))
    return {'x': x, 'z': z, 'y': y,
            'corrupt_prob': corrupt_prob,
            'cost': cost}


# %%


def test_flow_cap(target):
    import tensorflow as tf
    import tensorflow.examples.tutorials.mnist.input_data as input_data
    import matplotlib.pyplot as plt

    # %%
    # load MNIST as before
    # mnist = input_data.read_data_sets('MNIST_data', one_hot=True)
    # mean_img = np.mean(mnist.train.images, axis=0)
    # ae = autoencoder(dimensions=[784, 256, 64])
    if target == 'flowcap':
        flow_cap_dir = '/Local-Scratch/CMPT-884-workspace/flowCAP-II/result/feat/file-cell-countAdj.csv'
        data, sample_number_list = read_csv(flow_cap_dir)
    elif target == 'IMPC':
        IMPC_dir = '/Local-Scratch/CMPT-884-workspace/IMPC/result/P1/Sanger_SPLEEN/feat/file-cell-countAdj.csv'
        data, sample_number_list = read_csv(IMPC_dir)
    else:
        raise ValueError("wrong target")
    shape = data.shape
    feature_num = shape[1]
    data_num = shape[0]
    learning_rate_list = [0.0001, 0.001, 0.01]
    plt.figure(figsize=(7, 5))
    for learning_rate in learning_rate_list:
        start_time = time.time()
        ae = autoencoder(dimensions=[feature_num, feature_num / 2, feature_num / 3, feature_num / 4, feature_num / 6,
                                     feature_num / 8])

        # %%
        # learning_rate = 0.001
        optimizer = tf.train.AdamOptimizer(learning_rate).minimize(ae['cost'])

        # %%
        # We create a session to use the graph
        with tf.Session() as sess:
            sess.run(tf.global_variables_initializer())

            iteration_number_plot = []
            loss_plot = []

            # %%
            # Fit all training data
            batch_size = 100
            iteration_number = 100

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
                    _, l = sess.run([optimizer, ae['cost']], feed_dict={ae['x']: batch_data, ae['corrupt_prob']: [1]})
                    # Display logs per step
                    # if i % display_step == 0 or i == 1:
                l, y = sess.run([ae['cost'], ae['y']], feed_dict={ae['x']: data, ae['corrupt_prob']: [1]})
                print('iteration %i: Minibatch Loss: %f' % (i, l))

                loss_plot.append(l)
                iteration_number_plot.append(i)

            dimension_reduction = (sess.run([ae['z']], feed_dict={ae['x']: data, ae['corrupt_prob']: [0]}))[0]
        plt.plot(iteration_number_plot, loss_plot, linewidth=3, label='LR={0}'.format(learning_rate))
        print("---sdae for {0} {1} seconds ---".format(target, time.time() - start_time))

    plt.xlabel('Iteration', fontsize=18)
    plt.ylabel('Loss', fontsize=18)
    plt.legend(loc="best", fontsize=15)
    plt.savefig('./draw_result/{1}-sdae-loss-{0}-lr.png'.format(iteration_number, target))

    # save_csv_data(dir='./dr_results/{0}-sdae-{1}.csv'.format(target, int(feature_num / 8)), data=dimension_reduction,
    #               sample_number_list=sample_number_list,
    #               name='sdae')
    # draw_loss_function(iteration_number_plot, loss_plot, './draw_result/{0}-sdae-loss-100'.format(target))

    # for epoch_i in range(n_epochs):
    #     for batch_i in range(mnist.train.num_examples // batch_size):
    #         batch_xs, _ = mnist.train.next_batch(batch_size)
    #         train = np.array([img - mean_img for img in batch_xs])
    #         sess.run(optimizer, feed_dict={
    #             ae['x']: train, ae['corrupt_prob']: [1.0]})
    #     print(epoch_i, sess.run(ae['cost'], feed_dict={
    #         ae['x']: train, ae['corrupt_prob']: [1.0]}))

    # %%
    # Plot example reconstructions
    # n_examples = 15
    # test_xs, _ = mnist.test.next_batch(n_examples)
    # test_xs_norm = np.array([img - mean_img for img in test_xs])
    # recon = sess.run(ae['y'], feed_dict={
    #     ae['x']: test_xs_norm, ae['corrupt_prob']: [0.0]})
    # fig, axs = plt.subplots(2, n_examples, figsize=(10, 2))
    # for example_i in range(n_examples):
    #     axs[0][example_i].imshow(
    #         np.reshape(test_xs[example_i, :], (28, 28)))
    #     axs[1][example_i].imshow(
    #         np.reshape([recon[example_i, :] + mean_img], (28, 28)))
    # fig.show()
    # plt.draw()
    # plt.waitforbuttonpress()


if __name__ == '__main__':
    # test_flow_cap(target='IMPC')
    test_flow_cap(target='flowcap')
