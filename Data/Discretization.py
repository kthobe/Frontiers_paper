# Program for the discretization of Bioplex data
# Kirsten Thobe --- 10.08.2018

# Main program starts below line

import csv
import numpy as np

# function for importing data from CSV file and returning them as matrix
def CSV2Matrix(filename):
    print(f"Reading in CSV file {filename}")
    with open(filename) as csvfile:
            data_raw = csv.reader(csvfile, delimiter=',')
            # remove headers
            line1 = next(data_raw)
            first = line1[0]                # first line with headers
            print(f"Species in the dataset are {line1[1:]}")
            maxlen = len(line1)             # number of components
            list_ = []                      # create matrix of data set
            tpoints =[]
            for line in data_raw:
                if(line[0] == first):
                    break;
                else:
                    list_.append(line[1:maxlen])
                    tpoints.append(line[0])
    return( np.asarray(list_).astype(np.float), line1[1:], np.asarray(tpoints))

# function for calculating the mean or median for two datasets and discretizing the data
def Real2Bool(dataset1,dataset2,method):
    combidata = np.concatenate((dataset1,dataset2),axis=0)      # combine both datasets
    data = combidata.transpose()
    maxlen = data.shape[1]                                      # size of combined datasets
    lendata = dataset1.shape[0]                                 # size of one datasets
    result_array = np.empty((0,maxlen), int)
    meanvalues = []
    for row in data:
        if(method=='mean'):                                     # calculating mean or median per species
            result = np.mean(row)
        else:
            result = np.median(row)

        meanvalues = np.append(meanvalues, result)
        mean_array = np.array([])
        for i in range(maxlen):                                 # discretization by threshold
            if(row[i]>=result):
                boolval=1
                mean_array = np.append(mean_array, boolval)
            else:
                boolval=0
                mean_array = np.append(mean_array, boolval)
        result_array = np.append(result_array, [mean_array], axis=0)
    dis_data1 = np.int_(result_array[:,0:lendata].transpose())
    dis_data2 = np.int_(result_array[:,lendata:maxlen].transpose())
    print(f"Discretized values of species by {method} are:")
    print(meanvalues)
    return(dis_data1,dis_data2)
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------

# MAIN PROGRAM

# enter name of datasets to discretize, must be in pairs
filename_treated = 'Bp1851Sora.csv'         # options: 'Bp1851Sora.csv', 'Bp1851_2Sora.csv', 'Bp1257Sora.csv', 'Bp1257_2Sora.csv'
filename_untreated = 'Bp1851DMSO.csv'      # options: 'Bp1851DMSO.csv', 'Bp1851_2DMSO.csv', 'Bp1257DMSO.csv', 'Bp1257_2DMSO.csv'

# select discretization method
dis_method = 'median'                         # 'mean' or 'median'

datamatrix, species, tpoints = CSV2Matrix(filename_treated)
datamatrix2, species, tpoints = CSV2Matrix(filename_untreated)

# Annotate discretized data with measurement times
booldata1, booldata2 = Real2Bool(datamatrix,datamatrix2,dis_method)
annotate_data = np.column_stack((tpoints,booldata1))
annotate_data2 = np.column_stack((tpoints,booldata2))

print(f"Saving discretized data in new files called Dis{filename_treated} and Dis{filename_untreated}.")
np.savetxt('Dis'+filename_treated, annotate_data.astype(str), fmt='%s', delimiter='\t' ,header='\t'.join(species))
np.savetxt('Dis'+filename_untreated, annotate_data2.astype(str), fmt='%s', delimiter='\t', header=','.join(species))
