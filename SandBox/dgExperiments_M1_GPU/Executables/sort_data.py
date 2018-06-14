
import matplotlib.pyplot as plt

import re
import sys

if( len( sys.argv ) < 3 ):
    print('Please use 2 command line arguments')
    sys.exit()

thread_num = sys.argv[1]
mkl_thread_num = sys.argv[2]

fp = open('data.txt', 'r')

contents = fp.readlines()

fp.close()

output_order = ['Comp Time', 'numerical flux', 'left flux', 'Prep Time', \
        'second dgemm calls', 'Total time', 'update of result', \
        'right flux', 'first dgemm calls', \
        'ff', 'second num flux update', 'flux_x1_q', 'copy left/right', \
        'ef', 'primitive']

myData = {}

myCounts = {}

threads = [1, 2, 4, 8]


for line in contents:
    
    it = re.finditer('([a-zA-Z_0-9\s/]+):\s+(\d+.\d+)(E)?(-?\d{3})?', line)

    while( True ):
        try:
            match = next( it )
            tup = match.groups()
            name = tup[0].strip()
            fl = tup[1]
            exp = 0

            if( name == "InitializeFields" ):
                continue

            if( tup[3] != None ):
                exp = tup[3]

            if name not in myData.keys():
                myData[name] = 0.0

            if name not in myCounts.keys():
                myCounts[name] = 0

            myCounts[name] += 1
            myData[name] += float( fl ) * pow(10,float( exp ))
        except StopIteration:
            break

for key in myData.keys():
    myData[key] /= myCounts[key]
    
fp = open('output.txt', 'a+')

fp.write('{},{}\t'.format(thread_num,mkl_thread_num))

if( len(output_order) != len(myData.keys()) ):
    print('error with not including all keys')

for i in range(0,len(output_order)):
    st = output_order[i]

    if st not in myData.keys():
        print('HUGE ERROR, KEY NOT FOUND')

    fp.write('{0:<1.5f}\t'.format(myData[st]))

fp.write('\n')

fp.close()

