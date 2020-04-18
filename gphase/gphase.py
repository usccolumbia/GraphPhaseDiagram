import sys
import matplotlib.pyplot as plt
import matplotlib.colors
import phase_module
import argparse
import matlab.engine
import os
from argparse import RawTextHelpFormatter
import numpy as np
import peakdetect
import warnings
import copy

if __name__ == "__main__":
    
    # ./bin/GPhase.py -i data/xrd.csv -c data/composition.sv 
    # ./bin/GPhase.py -i xrd.csv -c composition.sv 
    # Call matlab function

    warnings.filterwarnings("ignore")
    parser = argparse.ArgumentParser(description="GPhase--An XRD Phase Mapping Program.\n University of South Carolina",
                formatter_class=RawTextHelpFormatter)
    parser.add_argument('-x', '--xrd', dest='xrd', help='x-ray diffraction csv file', default='../data/FePdGa_XRD.csv')
    parser.add_argument('-c', '--comp', dest='comp', help='composition csv file', default='../data/FePdGa_composition.csv')
    parser.add_argument('-k', "--K", dest='K', help='threshold for merging phases (default:1.5)', type=float, default=1.5)
    parser.add_argument('-b', '--background', dest='if_background_subtraction', help='background subtraction', type=int, default=0)
    # parser.add_argument('-v', "--version")
    args = parser.parse_args()

    # execute MATLAB peak finder code
    if args.if_background_subtraction == 0:
        eng = matlab.engine.start_matlab()
        bin_path = os.path.dirname(os.path.realpath(__file__))
        eng.cd(bin_path, nargout=0)
        print "finding peaks... takes a few minutes..."
        xrd_peak_path = eng.peakfinder(args.xrd, nargout=1)
    else:
        xrd_peak_path = args.xrd
        print "calculating..."

    # read data
    [xrd, two_theta, composition, label] = phase_module.read_data(xrd_peak_path, args.comp)
    sample_number = len(xrd)
    feature_number = len(xrd[0])

    # plt plot the samples
    # x = []
    # y = []
    # for i in range(sample_number):
    #     x.append(composition[i][0])
    #     y.append(composition[i][1])
    # plt.scatter(x, y, c=label, cmap=plt.cm.jet)
    # plt.show()

    # Resize for ternary figure
    # for row in composition:
    #     row[0] = row[0] * 5 / 3
    #     row[1] = row[1] * 5 / 3
    #     row[2] = (row[2] - 40) * 5 / 3

    # Construct neighbors
    [neighbor_list_original, neighbor_list, original_comp, coordinate, x_coordinate, y_coordinate] = \
    phase_module.construct_neighbor(composition)
    # four_nearest_neighbor = phase_module.nearest_neighbor(composition, 4)
    # eight_nearest_neighbor = phase_module.nearest_neighbor(composition, 8)

    # Plot Delaunay triangulation
    position = None
    scale = 100.0
    font_size = 20
    text_content = ["(at.%)", "(at.%)", "(at.%)"]
    text_position = [(0,-7.5,107.5), (102.5,-7.5,5), (-5,102,3)]
    ternary_data = [scale, position, font_size, text_content]
    # for data in [neighbor_list_original, neighbor_list]:
    #     [figure, tax] = phase_module.ternary_figure(ternary_data)
    #     tax.scatter(original_comp, marker='o', edgecolor='w', s=40, color='red')
    #     for row in data:
    #         plt.plot([coordinate[row[0]][0], coordinate[row[1]][0]], [coordinate[row[0]][1], coordinate[row[1]][1]], 'b')
    #         tax.line(p1=(composition[row[0]][0],composition[row[0]][1],composition[row[0]][2]), p2=(composition[row[1]][0],composition[row[1]][1],composition[row[1]][2]))
    #     plt.show()
    #     tax.close()

    if args.if_background_subtraction != 0:
        # calculate deviation
        std = np.std(xrd.tolist(), axis=0)
        std = std.tolist()
        std = phase_module.back_sub(std, neighbor=2, threshold=0.5, fitting_degree=50, if_plot=0, two_theta=two_theta)
        std = phase_module.back_sub(std, neighbor=2, threshold=0.8, fitting_degree=50, if_plot=0)
        for i in range(feature_number):
            if std[i] < 0:
                std[i] = 0

        # use python code peakdetect to detect peaks
        _max, _min = peakdetect.peakdetect(std, range(feature_number), lookahead=3, delta=0.35)
        # _max, _min = peakdetect.peakdetect_zero_crossing(std[1:100], window=5)
        max_list = []
        for i in range(len(_max)):
            max_list.append(_max[i][0])
        min_list = []
        for i in range(len(_min)):
            min_list.append(_min[i][0])
        xm = [p[0] for p in _max]
        ym = [p[1] for p in _max]
        xn = [p[0] for p in _min]
        yn = [p[1] for p in _min]

        # find peak area for deviation
        if min_list[0] < max_list[0]:
            shift = 1
        else:
            shift = 0
        filter_list = [0] * feature_number
        for i in range(len(max_list)):
            peak_position = max_list[i]
            filter_area = [peak_position, peak_position]
            if_edge = 0
            while if_edge == 0:
                filter_area[0] -= 1
                filter_area[1] += 1
                if std[filter_area[0]] <= 0 or std[filter_area[1]] <= 0:
                    if_edge = 1
                if i == 0 and shift == 0:
                    if filter_area[1] >= min_list[i + shift]:
                        if_edge = 1
                elif i == len(max_list) - 1 and min_list[-1] < max_list[-1]:
                    if filter_area[0] <= min_list[i - 1 + shift]:
                        if_edge = 1
                else:
                    if filter_area[0] <= min_list[i - 1 + shift] or filter_area[1] >= min_list[i + shift]:
                        if_edge = 1
            # Minimum peak area size
            if filter_area[1] - filter_area[0] >= 8:
                filter_list[filter_area[0]:(filter_area[1] + 1)] = [1] * (filter_area[1] - filter_area[0] + 1)
        std_peak = [a * b for a, b in zip(std, filter_list)]

        # Plot std figure before and after peak filter
        # plt.subplot(211)
        # plt.plot(two_theta, std)
        # plt.xlim(two_theta[0], two_theta[-1])
        # for i in range(len(xm)):
        #     xm[i] = two_theta[xm[i]]
        # for i in range(len(xn)):
        #     xn[i] = two_theta[xn[i]]
        # plt.plot(xm, ym, 'r*')
        # plt.plot(xn, yn, 'go')
        # plt.subplot(212)
        # plt.plot(two_theta, std_peak)
        # plt.xlim(two_theta[0], two_theta[-1])
        # plt.show()
        # plt.close()

        # Plot std figure after filter and filter figure
        # plt.figure(dpi=150)
        # plt.xlabel(r'$2\theta\,Angle (Deg.)$', fontsize=20)
        # plt.ylabel(r'$Intensity$', fontsize=20)
        # plt.plot(two_theta, std_peak)
        # plt.xlim(two_theta[0], two_theta[-1])
        # plt.savefig('5-3.png', format='png')
        # plt.close()
        # plt.figure(dpi=150)
        # plt.xlabel(r'$2\theta\,Angle (Deg.)$', fontsize=20)
        # plt.ylabel(r'$Intensity$', fontsize=20)
        # plt.plot(two_theta, filter_list)
        # plt.xlim(two_theta[0], two_theta[-1])
        # plt.ylim(-0.2, 2.0)
        # plt.savefig('5-4.png', format='png')
        # plt.close()

        # apply filter to each samples
        xrd_filterdata = []
        xrd_temp = copy.deepcopy(xrd.tolist())
        if_plot = 0
        for i in range(sample_number):
            xrd_temp[i] = phase_module.back_sub(xrd_temp[i], neighbor=2, threshold=0.5, fitting_degree=50, if_plot=0,
                                         two_theta=two_theta)
            xrd_temp[i] = phase_module.back_sub(xrd_temp[i], neighbor=2, threshold=0.8, fitting_degree=50, if_plot=0)
            for j in range(feature_number):
                if xrd_temp[i][j] < 0:
                    xrd_temp[i][j] = 0
            xrd_filterdata.append([a * b for a, b in zip(xrd_temp[i], filter_list)])
        # xrd_curve = []
        for sample in range(sample_number):
            for feature in range(feature_number):
                # xrd_curve.append(xrd_filterdata[sample][feature])
                if xrd_filterdata[sample][feature] > 10:
                    xrd_filterdata[sample][feature] = 1
                else:
                    xrd_filterdata[sample][feature] = 0
        xrd = xrd_filterdata
    else:
        # normalization
        for sample in range(sample_number):
            max_data = max(xrd[sample])
            if max_data != 0:
                xrd[sample][:] = [x/max_data for x in xrd[sample]]


    neighbor_list = phase_module.euclid_distance(data=xrd, neighbor_list=neighbor_list)

    # edge weight normalization
    max_neighbor_list = [max(b) for b in zip(*neighbor_list)][2]
    for row in range(len(neighbor_list)):
        neighbor_list[row][2] /= max_neighbor_list

    # output all neighbors and their distances
    # with open('neighbor.csv', 'wb') as csvfile:
    #     spamwriter = csv.writer(csvfile)
    #     for row in neighbor_list:
    #         spamwriter.writerow((row[0], row[1], row[2]))

    # image segmentation
    num_nodes = sample_number - 1
    min_size = 1
    [forest, prediction] = phase_module.graph_based_segmentation(neighbor_list, num_nodes, args.K, min_size, phase_module.threshold)
    print 'Number of components: %d' % forest.num_sets

    # Evaluate results
    phase_module.result_evaluation(label, prediction)

    # Plot results
    len_truth = max(label)
    len_predicted = max(prediction)
    [figure, tax] = phase_module.ternary_figure(ternary_data)
    tax.scatter(original_comp, marker='o', c=label, s=30, norm=matplotlib.colors.LogNorm(), cmap=plt.cm.jet)
    tax.savefig("../figure/" + "labeled.png", format="png")
    # tax.show()

    [figure, tax] = phase_module.ternary_figure(ternary_data)
    tax.scatter(original_comp, marker='o', c=prediction, s=30, norm=matplotlib.colors.LogNorm(), cmap=plt.cm.jet)
    tax.savefig("../figure/" + "prediction.png", format="png")
    # tax.show()

    print "Calculation finished."