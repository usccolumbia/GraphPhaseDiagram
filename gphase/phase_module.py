"""
Data analysis modules
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import csv
import sys
import ternary
import copy
from math import acos, pi
from graph import segment_graph, build_graph
from random import random
from scipy.spatial import Delaunay

def threshold(size, const):
    return const / size

def read_data(xrd, comp):
    """
    read composition data and x-ray diffraction data from csv file.

    Parameters
    ----------

    """
    # read composition file and XRD file
    xrd_data = np.genfromtxt(xrd, delimiter=',')
    comp_data = np.genfromtxt(comp, delimiter=',')
    two_theta = xrd_data[0, :]
    xrd_data = xrd_data[1:, :]
    label = comp_data[:, 3]
    comp_data = comp_data[:, 0:3]
    return [xrd_data, two_theta, comp_data, label]


def back_sub(data, neighbor, threshold, fitting_degree, if_plot, two_theta=None):
    """
    Use polynomial fitting to do background subtraction for curves.

    """
    feature_number = len(data)
    x = range(feature_number)
    y = copy.deepcopy(data)

    # Calculate difference between spectrum points
    y_diff = [0] * feature_number

    # Calculate neighbor difference
    for i in range(neighbor, feature_number-neighbor):
        for j in range(1, neighbor+1):
            y_diff[i] = y_diff[i] + abs(y[i] - y[i+j]) + abs(y[i] - y[i-j])
    for i in range(neighbor):
        y_diff[i] = float("Inf")
        y_diff[feature_number-i-1] = float("Inf")

    # Calculate subtraction threshold, delete points which past the threshold
    y_diff_sort = sorted(y_diff)
    min_threshold = y_diff_sort[int(round(len(y_diff_sort)*threshold))-1]
    for i in range(feature_number-1, -1, -1):
        if y_diff[i] > min_threshold:
            del x[i]
            del y[i]

    if if_plot != 0:
        x_twotheta = []
        for i in range(len(x)):
            x_twotheta.append(two_theta[x[i]])
        x_twotheta = np.array(x_twotheta)

    if if_plot != 0:
        # Original figure
        plt.figure(dpi=150)
        plt.xlabel(r'$2\theta\,Angle (Deg.)$', fontsize=20)
        plt.ylabel(r'$Intensity$', fontsize=20)
        plt.plot(two_theta, data)
        plt.xlim(two_theta[0], two_theta[-1])
        # plt.savefig('', format='png')
        # diff figure with scattered interest point
        plt.figure(dpi=150)
        plt.xlabel(r'$2\theta\,Angle (Deg.)$', fontsize=20)
        plt.ylabel(r'$Diff\,Intensity$', fontsize=20)
        plt.plot(two_theta, y_diff)
        plt.xlim(two_theta[0], two_theta[-1])
        for i in range(len(x)):
            plt.scatter(two_theta[x[i]], y_diff[x[i]], s=15, marker='.', color='red')
        # plt.savefig('g', format='png')
        # Original figure with scattered interest point
        plt.figure(dpi=150)
        plt.xlabel(r'$2\theta\,Angle (Deg.)$', fontsize=20)
        plt.ylabel(r'$Intensity$', fontsize=20)
        plt.plot(two_theta, data)
        plt.xlim(two_theta[0], two_theta[-1])
        plt.scatter(x_twotheta, y, s=15, marker='.', color='red')
        # plt.savefig('', format='png')

    # Polynomial fit, p will be the polynomial function
    x = np.array(x)
    y = np.array(y)
    z = np.polyfit(x, y, fitting_degree)
    p = np.poly1d(z)

    # backup data for plot
    data_backup = copy.deepcopy(data)

    # subtract data with polynomial curve
    for i in range(neighbor, feature_number-neighbor):
        data[i] = data[i] - p(i)
    for i in range(neighbor):
        data[i] = p(i)
        data[feature_number-i-1] = p(feature_number-i-1)
        # data[i] = 0
        # data[feature_number-i-1] = 0

    # plot or not
    if if_plot != 0:
        #plt.scatter(x,y)
        plt.figure(dpi=150)
        plt.xlabel(r'$2\theta\,Angle (Deg.)$', fontsize=20)
        plt.ylabel(r'$Intensity$', fontsize=20)
        x1 = 0 + 20
        x2 = feature_number - 20
        plt.plot(two_theta[x1:x2], data_backup[x1:x2], color='blue', linewidth=1.5)
        plt.plot(x_twotheta, p(x), color='red', linewidth=1.5)
        plt.plot(two_theta[x1:x2], data[x1:x2], color='green', linewidth=1.5)
        plt.xlim(two_theta[x1], two_theta[x2-1])
        # plt.savefig('', format='png')
        plt.close()

    return data


def construct_neighbor(comp_data):
    """
	Use Delaunary Triagulation to construct neighbor list.

	"""
    # Convert original composition to coordinate points data
    original_comp = []
    for point in comp_data:
        x = point[0]
        y = point[1]
        z = point[2]
        original_comp.append((x,y,z))
    xs, ys = ternary.helpers.project_sequence(original_comp, permutation=None)
    coordinate = []
    for num in range(len(xs)):
        coordinate.append([xs[num], ys[num]])
    coordinate = [[float(column) for column in row] for row in coordinate]
    coordinate = np.array(coordinate)

    # Original triangulation
    tri = Delaunay(coordinate)
    original_neighbors = []
    for row in tri.simplices:
        original_neighbors.append(sorted([row[1], row[2]]))
        original_neighbors.append(sorted([row[0], row[2]]))
        original_neighbors.append(sorted([row[0], row[1]]))

    # Remove duplicate neighbors
    neighbor_list = []
    for row in original_neighbors:
        if row not in neighbor_list:
            neighbor_list.append(row)

    neighbor_list_original = neighbor_list[:]

    # Modify triangulation to construct the new neighbor_list
    for row in tri.simplices:
        side1 = ((coordinate[row[0]][0] - coordinate[row[1]][0]) ** 2 +  (coordinate[row[0]][1] - coordinate[row[1]][1]) ** 2) ** 0.5
        side2 = ((coordinate[row[0]][0] - coordinate[row[2]][0]) ** 2 +  (coordinate[row[0]][1] - coordinate[row[2]][1]) ** 2) ** 0.5
        side3 = ((coordinate[row[1]][0] - coordinate[row[2]][0]) ** 2 +  (coordinate[row[1]][1] - coordinate[row[2]][1]) ** 2) ** 0.5
        if acos((side1**2 + side2**2 - side3**2)/(2*side1*side2)) > 0.7 * pi:
            neighbor_list.remove(sorted([row[1], row[2]]))
        if acos((side1**2 + side3**2 - side2**2)/(2*side1*side3)) > 0.7 * pi:
            neighbor_list.remove(sorted([row[0], row[2]]))
        if acos((side2**2 + side3**2 - side1**2)/(2*side2*side3)) > 0.7 * pi:
            neighbor_list.remove(sorted([row[0], row[1]]))

    return [neighbor_list_original, neighbor_list, original_comp, coordinate, xs, ys]


def ternary_figure(ternary_data):
    """
    Construct ternay triagulation figure.

    """
    [scale, position, font_size, text_content] = ternary_data
    figure, tax = ternary.figure(scale=scale)
    tax.boundary(linewidth=1.5)
    tax.gridlines(multiple=20, color="blue")
    tax.left_axis_label(text_content[0], fontsize=font_size, offset = 0.12)
    tax.right_axis_label(text_content[1], fontsize=font_size, offset = 0.12)
    tax.bottom_axis_label(text_content[2], fontsize=font_size, offset = -0.02)
    tax.ticks(axis='l', ticks=["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"], offset=0.022)
    tax.ticks(axis='b', ticks=["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"], offset=0.022)
    tax.ticks(axis='r', ticks=["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"], offset=0.022)
    tax.clear_matplotlib_ticks()
    return [figure, tax]
   
def euclid_distance(data, neighbor_list):
    """
    calculate Euclidean Distances between neighbors, add to neighbor_list

    """
    for line in range(len(neighbor_list)):
        i = neighbor_list[line][0]
        j = neighbor_list[line][1]
        summary = 0.0
        max_sample1 = max(data[i])
        max_sample2 = max(data[j])
        max_sample = max(max_sample1, max_sample2)
        if max(max_sample1, max_sample2) == 0:
            neighbor_list[line].append(0)
        else:
            for feature in range(len(data[0])):
                # if peakcenter_location[i][feature] - peakcenter_location[j][feature] > 5:
                summary += (data[i][feature]/max_sample - data[j][feature]/max_sample)** 2
            neighbor_list[line].append(summary ** 0.5)
    return neighbor_list

def graph_based_segmentation(neighbor_list, num_nodes, K, min_size, threshold):
    xcoor = []
    ycoor = []
    diff = []
    length = 0
    for i in range(len(neighbor_list)):
        xcoor.append(neighbor_list[i][0])
        ycoor.append(neighbor_list[i][1])
        diff.append(neighbor_list[i][2])
    graph = build_graph(length, diff, xcoor, ycoor)
    forest = segment_graph(graph, num_nodes + 1, K, min_size, threshold)

    # output prediction results
    prediction = []
    for x in xrange(num_nodes + 1):
        prediction.append(forest.find(x))

    # prediction clusters rename
    name_set = []
    for name in prediction:
        if name not in name_set:
            name_set.append(name)
    rename = {}
    for i in range(len(name_set)):
        rename[name_set[i]] = i + 1
    for i in range(len(prediction)):
        prediction[i] = rename[prediction[i]]

    # Write results to files and save
    # with open('prediction','w') as f:
    #     for i in prediction:
    #         f.write(str(i) + '\n')

    return [forest, prediction]

def result_evaluation(label, prediction):
    """
    Cluster result evaluation.

    """
    tp = 0.0
    tn = 0.0
    fp = 0.0
    fn = 0.0
    truth = 0.0
    classify = 0.0
    for i in range(len(label)-1):
        for j in range(i + 1, len(label)):
            if label[i] == label[j] and prediction[i] == prediction[j]:
                tp += 1
            if label[i] != label[j] and prediction[i] != prediction[j]:
                tn += 1
            if label[i] == label[j] and prediction[i] != prediction[j]:
                fn += 1
            if label[i] != label[j] and prediction[i] == prediction[j]:
                fp += 1
            if label[i] == label[j]:
                truth += 1
            if prediction[i] == prediction[j]:
                classify += 1
    # Calculate precision and recall
    print "precision = %f" % (tp / (tp + fp))
    print "recall = %f" % (tp / (tp + fn))
    print "accuracy = %f" % ((tp + tn) / (tp + tn + fp + fn))
    print "mcc = %f" % ((tp * tn - fp * fn) / (((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5))

def nearest_neighbor(comp_data, neighbor_num):
    """"
    Nearest neighbors

    """
    neighbor_list = []
    original_comp = []
    for point in comp_data:
        x = point[0]
        y = point[1]
        z = point[2]
        original_comp.append((x,y,z))
    xs, ys = ternary.helpers.project_sequence(original_comp, permutation=None)

    for i in range(len(comp_data)):
        weight = []
        for j in range(len(comp_data)):
            if i != j:
                distance = ((xs[i] - xs[j]) ** 2 + (ys[i] - ys[j]) ** 2)**0.5
                weight.append([i, j, distance])
        weight = sorted(weight, key = lambda x: x[2])
        for row in range(neighbor_num):
            neighbor_list.append(weight[row])
    neighbor_list = sorted(neighbor_list, key = lambda x: x[2])

    duplicate_neighbor_list = neighbor_list
    neighbor_list = []
    for row in range(len(duplicate_neighbor_list) - 1):
        if duplicate_neighbor_list[row][2] != duplicate_neighbor_list[row + 1][2]:
            neighbor_list.append(duplicate_neighbor_list[row])
    neighbor_list.append(duplicate_neighbor_list[-1])

    return neighbor_list