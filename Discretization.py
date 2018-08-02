# Discretization of Bioplex data

import csv

with open('Bp1257DMSO.csv') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        # remove headers
        #Bp1257DMSO = data[]
        x=[]
        for row in data:
            x.append(row[1])
            print(x)

# mean_array = np.array([])
# median_array = np.array([])
#
# for line in Bp1851exp1:
#     result = do_stuff(line)
#     result_array = np.append(result_array, result)
