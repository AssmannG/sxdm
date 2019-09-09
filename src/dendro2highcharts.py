#!/usr/bin/env python

# Usage in the code

# from dendro2highcharts import dendro2highcharts
# scipy_dendro -> is output of scipy.cluster.hierarchy
# dendro_formatted = dendro2highcharts(scipy_dendro)

import json

def dendro2highcharts(dendro_dict):
    '''Function that converts output of scipy.cluster.hierarchy.dendrogram
       to data that can be easily plot using Highcharts.js scatterplot

       => Input:
            dendro_dict: dictonary outputed by scipy.cluster.hierarchy.dendrogram
    '''
    try:
        X = dendro_dict['icoord']
        Y = dendro_dict['dcoord']
        leaves = dendro_dict['ivl']
        colors = dendro_dict['color_list']
    except KeyError as e:
        raise Exception('dendro2highcharts: Required key not found in dendro_dict: {}'.format(e))
    sorted_by_colors = {}
    leaves_indices = []
    x_labels = {}

    # Organize a dictionary with [x,y] points assigned to each respective color.

    for i in range(len(X)):
        x = X[i]
        y = Y[i]
        color = colors[i]
        if color not in sorted_by_colors.keys():
            sorted_by_colors[color] = []
        for j in range(4):
            elem = [x[j], y[j]]
            sorted_by_colors[color].append(elem)
            if j == 3:
                sorted_by_colors[color].append(None)
            if y[j] == 0:
                leaves_indices.append(x[j])

    # Sort the labels in the ascending orders
    leaves_indices.sort()

    series = []
    # Convert colors sort to highcharts series
    for k in sorted_by_colors:
        d = sorted_by_colors[k]
        series.append({'data':d})

    # Associates leaf labesl with its position on x_axis
    for i, x in enumerate(leaves_indices):
        l = leaves[i]
        x = '{}{}'.format('x',int(x))
        x_labels[x] = l


    output = {}
    output['series'] = series
    output['labels'] = x_labels

    #output_json = json.dumps(output)

    return output


if __name__ == '__main__':

    # Load hcluster matrix from hclust from the file and genarate dendrogram
    import scipy.cluster.hierarchy as sch
    import hclust
    import numpy as np
    hc = hclust.h
    hca = np.array(hc)
    d = sch.dendrogram(hca, p=10)

    highcharts_data = dendro2highcharts(d)

    print (highcharts_data)
