#!usr/bin/env python3

import logging
import argparse
import pandas as pd
pd.set_option('display.max_columns', 500)
import glob
from sklearn.cluster import KMeans
log = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description='Find minE file from a trafl-file')
    parser.add_argument ('-i', '--input', help='Path to Inputfiles')
    parser.add_argument ('-n', '--number', type=int, help='Number of "--save-n-best" in ernwin')
    parser.add_argument ('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument ('-c', '--cluster', type=int, help='number of clusters')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    dict_elements=dict()
    df_elements = pd.DataFrame()
    df_indexelements=pd.DataFrame()
    elementindex= list()
    dict_elementindex=dict()
    column_names=list()  #column names for dataframe
    index_names=list() #index names for dataframe

    for n in range(int(args.number)):
        with open (args.input+'best'+str(n)+'.coord', 'r') as FILE:
            key='best'+str(n) #best ernwin structure file
            column_names.append(key)

            for ln in FILE:
                if ln.startswith('sampled'): #indicate a new element in a best.coord file
                    samples=(ln[8:].split())

                    if samples[0] not in index_names:
                        index_names.append(samples[0])

                    if samples[1] not in elementindex: #sample[1] e.g. 4YCP_B:m_3
                        elementindex.append(samples[1])

                    value_elements = ','.join(samples[1::1]) # 4YCP_B:m_3,2,4,5
                    value_elementsindex = elementindex.index(samples[1])
                    #samples[0] e.g. s0,s1,m0,i0
                    if key in dict_elements: #add element to existing best.coord
                        dict_elements[key].append([samples[0],value_elements])
                        dict_elementindex[key].append([samples[0],value_elementsindex])
                    else: #new best.coord and add new element
                        dict_elements[key]=[[samples[0],value_elements]]
                        dict_elementindex[key]=[[samples[0],value_elementsindex]]

                else:
                    continue

    #create empty dataframes
    df_elements = pd.DataFrame(columns=column_names, index=index_names)
    df_indexelements = pd.DataFrame(1000,columns=index_names, index=column_names)

    for k, value in  dict_elements.items():
        for v in value:
            df_elements.at[v[0],[k]] = v[1]

    for k, value in dict_elementindex.items():
        for v in value:
            df_indexelements.at[[k],v[0]] = v[1]

    log.debug(df_elements)
    log.debug(df_indexelements)

    kmeans = KMeans(n_clusters=int(args.cluster))
    y = kmeans.fit_predict(df_indexelements)
    df_indexelements['cluster'] = y

    df_indexelements.index.name = 'Ernwinstructure'

    df_indexelements.to_csv(args.input+'clusterdataframe.csv', header=True)
    df_elements.to_csv(args.input+'usedelements.csv',header=True)

    grouped_df = df_indexelements.groupby('cluster')

    log.debug('grouped dataframes')
    log.debug(grouped_df)

    ernwinstructurelist = list()
    count = int(0)

    while count < args.cluster:
        log.debug('this is the count {}'.format(count))
        for cluster, df_cluster in df_indexelements.groupby('cluster'):
            log.debug('bestcluster')
            log.debug(df_cluster)
            log.debug('beststructure')
            best = (df_cluster.head(1))
            log.debug(best)

            number = best.index[0]
            log.debug('number')
            log.debug(number)
            ernwinstructurelist.append(number)

            outputfile = args.input+'cluster'+str(count)+'.csv'
            df_cluster.to_csv(outputfile, columns=[], header=False)

            count += 1

    print(ernwinstructurelist)

if __name__ == "__main__":
    main()
