import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

label_size = 15
mpl.style.use('seaborn')
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size


def save_csv_data(dir, data, sample_number_list, name):
    with open(dir, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        feature_name_list = []
        for i in range(0, len(data[0])):
            feature_name_list.append('{0}-{1}'.format(name, i + 1))
        feature_name_list.insert(0, 'sample-name')
        writer.writerow(feature_name_list)

        for index in range(0, len(data)):
            row = data[index].tolist()
            row.insert(0, sample_number_list[index])
            # np.insert(row, 0, sample_number_list[index])
            writer.writerow(row)


def read_csv(dir):
    data_array_all = []
    sample_number_list = []
    jump_row_flag = True
    with open(dir, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if jump_row_flag:
                jump_row_flag = False
                continue
            feature = map(float, row[1:])
            sample_number_list.append(row[0])
            data_array_all.append(np.asarray(feature))
    return np.asarray(data_array_all), sample_number_list


def draw_loss_function(x, y, save_name):
    plt.figure(figsize=(6, 5))
    plt.plot(x, y, linewidth=3)
    plt.xlabel('Iteration', fontsize=18)
    plt.ylabel('Loss', fontsize=18)
    plt.savefig(save_name)
    # plt.show()


