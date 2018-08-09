# Discretization of Bioplex data

import csv
import numpy as np

def CSV2Matrix(filename):
    with open(filename) as csvfile:
            data_raw = csv.reader(csvfile, delimiter=',')
            # remove headers
            line1 = next(data_raw)
            first = line1[0]                # first line with headers
            maxlen = len(line1)             # number of components
            list_ = []                      # create matrix of data set
            for line in data_raw:
                if(line[0] == first):
                    break;
                else:
    #            while(line[0] == first):
                    list_.append(line[1:maxlen])
    return(list_)


def Real2Bool(data,data2):
    print(data)
    print(data2)
    combidata = np.concatenate([data],[data2])
    print(combidata)
    maxlen = data.shape[1]
    result_array = np.empty((0,maxlen), int)
    result_array2 = np.empty((0,maxlen), int)
    #median_array = np.array([])
    for row in data:
        result = np.mean(row)
        mean_array = np.array([])
        for i in range(maxlen):
            if(row[i]>=result):
                boolval=1
                mean_array = np.append(mean_array, boolval)
            else:
                boolval=0
                mean_array = np.append(mean_array, boolval)
        print(mean_array)
        result_array = np.append(result_array, [mean_array], axis=0)
    return(result_array)

datamatrix = np.asarray(CSV2Matrix('Bp1851Sora.csv')).astype(np.float).astype(np.float)
datamatrix2 = np.asarray(CSV2Matrix('Bp1851DMSO.csv')).astype(np.float).astype(np.float)

combidata = np.append([datamatrix],[datamatrix2])
print(combidata)

#booldata = Real2Bool(datamatrix.transpose(),datamatrix2.transpose())
#print(booldata.transpose())
