# Discretization of Bioplex data

import csv
import numpy as np


with open('Bp1257DMSO.csv') as csvfile:
        data_raw = csv.reader(csvfile, delimiter=',')
        # remove headers
        data = data_raw[2-6,2-12]
        print(data)
        #Bp1257DMSO = data[]
        # x=[]
        # for row in data:
        #     x.append(row[1])
        #     print(x)

        # mean_array = np.array([])
        # median_array = np.array([])
        #
        # for row in data:
        #     print(row)
        #     result = np.mean(row)
            # mean_array = np.append(result_array, result)
            # print(result_array)
