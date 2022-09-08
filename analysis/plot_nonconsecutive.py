#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings
import glob, os
import io
import seaborn as sns
import matplotlib.ticker as plticker
warnings.simplefilter("ignore")
from matplotlib.legend import _get_legend_handles_labels

#Variable Settings
RNA = "1zciAB" #BASENAME
direct= "tests/" #PATH
clusters=10 #CLUSTER
exsteps=20 #ROUNDS
expansionsteps = list(range(3,20,2)) #start, end, interval
consecutive_run = "target" #for outputname: consecutive, nonconsecutive, target
consecutive = "consecutive"  #for outputname: consecutive, nonconsecutive,
count_type = 'nr' #choose between 'nr' or 'percent'
count_label = 'structure count [nr]'
energy_type = 'energy_min' #choose between:'energy_median', 'energy_mean', 'energy_min'
constraint_type = 'dif_median' #choose between:'dif_median', 'dif_mean', 'dif_min'

multipleruns = [direct+"runtime1",
                direct+"runtime2"]
                #if you want to plot more then one runtime at once

stepsize1 = ["5000", "10-000"] #sample runtimes

stepsize2 = ["expand_long00", "expand"] #sample runtimes
nrs = 5 #SIMROUND
seeds = 5 #nr. of SimRNA seeds


#colouring
blue = sns.light_palette("darkblue", as_cmap=True)
blue.set_under(color="w", alpha=None)

green = sns.light_palette("darkgreen", as_cmap=True)
green.set_under(color="w", alpha=None)

red = sns.light_palette("firebrick", reverse=True, as_cmap=True)
red.set_over(color="w", alpha=None)

brown = plt.cm.get_cmap('pink_r')


filename = list()
fclust=list()
df_sum = list()
subplotsname= list()
controll=list()

for p, path in enumerate(multipleruns):

    controll.append(str(consecutive+stepsize1[p]))
    #df_sum.append(str("df"+str(consecutive+stepsize[p])

    for runs in range(exsteps): # could be also lower/higher - if for one is no expansion possible
        run= str(runs)

        for cluster in range(clusters): #collect for each expansion step/cluster
            cl=str(cluster)

            for nr in range(nrs): #5  pdbs as starting points
                n=str(nr+1)

                for seed in range (seeds): #5 different seeds per pdb
                    #each single file
                    seed= str(seed+1)
                    file = str(path+"/cluster"+cl+"/"+run+"/"+
                               n+"/"+RNA+'c'+cl+"_"+stepsize2[p]+"_"+run+"_"+seed+".csv")

                    try:
                        df= pd.read_csv(file, sep="\t",skiprows=[i for i in range(1,103)])
                        df['stepsize'] = str(stepsize1[p])
                        df['constraint'] = int(expansionsteps[runs])
                        df['expansionstep'] = int(runs)
                        df['run'] = int(nr+1)
                        df['cluster'] = int(cluster)
                        df['seed'] = seed
                        fclust.append(df)
                        controll.append([file,int(stepsize1[p]),int(expansionsteps[runs]),
                                                  int(cluster),int(runs),int(nr+1),seed])

                    except:
                        #print("Probably there were no further expansion possible - please check")
                        #print(file)
                        pass

    df_c = pd.concat(fclust)


    df_sum.append(df_c)

df_consecutive10 = pd.concat(df_sum)

df_consecutive10["energy_diff"]= df_consecutive10["energy_values_plus_restraint_score"] \
                                        - df_consecutive10["energy_value"]

group1_collect = list()
group2_collect = list()



#plot all in one
for step in stepsize1:

    step_df= df_consecutive10[(df_consecutive10["stepsize"] == step)]

    title1 = (RNA+"_"+str(clusters)+consecutive+"_"+consecutive_run+"_"+str(step)+"_"+count_type)
    title2 = (RNA+"_"+str(clusters)+consecutive+"_"+consecutive_run+"_consecdif_"+str(step)+count_type)

    for run in range(exsteps):
        new_df = step_df[(step_df["constraint"] == expansionsteps[run])]
        group1 = new_df.groupby(["interaction_countbp"]).agg({'interaction_countbp': ['count'],
                                                            'energy_value': ['min', 'median', 'mean'],
                                                            'energy_diff':['min','max','median','mean']}
                                                    ).reset_index()
        group1.columns =['interaction_bp', 'nr',
                          'energy_min', 'energy_median','energy_mean',
                          'dif_min','dif_max', 'dif_median','dif_mean']
        group1['percent'] = (group1.nr /group1.nr.sum())*100
        group1['constraint'] = expansionsteps[run]
        group1.sort_values(by=['interaction_bp'])
        group1 = group1.astype({'energy_min': float, 'energy_median': float,'energy_mean': float})
        #group1 = group1.round({'energy_min': 0, 'energy_median': 0,'energy_mean': 0})
        group1_collect.append(group1)

        group2 = new_df.groupby(["len_interaction"]).agg({'len_interaction': ['count'],
                                                            'energy_value': ['min', 'median', 'mean'],
                                                            'energy_diff':['min','max','median','mean']}
                                                         ).reset_index()
        group2.columns =['interaction_length', 'nr',
                          'energy_min', 'energy_median','energy_mean',
                          'dif_min','dif_max', 'dif_median','dif_mean']
        group2['percent'] = (group2.nr /group2.nr.sum())*100
        group2['constraint'] = expansionsteps[run]
        group2.sort_values(by=['interaction_length'])
        group2 = group2.astype({'energy_min': float, 'energy_median': float,'energy_mean': float})
        group2_collect.append(group2)

    heatmapdata1 = pd.concat(group1_collect)
    heatmapdata2 = pd.concat(group2_collect)

    #CONSECUTIVE
    heat_intl0 = pd.pivot_table(heatmapdata1,
                            index='interaction_bp',
                            columns='constraint',
                            values= count_type,
                            aggfunc=np.sum,
                            fill_value=0)

    print(step)
    print(heat_intl0)
    heat_intl1 = np.dstack((heat_intl0,
                            np.zeros_like(heat_intl0),
                            np.zeros_like(heat_intl0))
                            ).reshape(heat_intl0.shape[0],-1)

    heat_intl12 = np.dstack((heat_intl0,
                            np.zeros_like(heat_intl0))
                            ).reshape(heat_intl0.shape[0],-1)


    heat_emin1 = pd.pivot_table(heatmapdata1,
                                 index='interaction_bp',
                                 columns='constraint',
                                 values= energy_type,
                                 aggfunc=np.sum,
                                 fill_value=0)

    heat_emin1 =np.dstack((np.zeros_like(heat_emin1),
                        heat_emin1,
                        np.zeros_like(heat_emin1))
                        ).reshape(heat_emin1.shape[0],-1)

    heat_edif1 = pd.pivot_table(heatmapdata1,
                               index='interaction_bp',
                               columns='constraint',
                               values= constraint_type,
                               aggfunc=np.sum,
                               fill_value=0)

    heat_edif1 =np.dstack((np.zeros_like(heat_edif1),
                        np.zeros_like(heat_edif1),
                        heat_edif1)
                        ).reshape(heat_edif1.shape[0],-1)

    print("constrain min {} constrain max {}".format(heatmapdata1[constraint_type].min()
                                                    ,heatmapdata1[constraint_type].max()))

    print("energy min {} energy max {}".format(heatmapdata1[energy_type].min()
                                                    ,heatmapdata1[energy_type].max()))

    #NON-CONSECUTIVE
    heat_intl2 = pd.pivot_table(heatmapdata2,
                            index='interaction_length',
                            columns='constraint',
                            values=count_type,
                            aggfunc=np.sum,
                            fill_value=0)

    heat_intl2 = np.dstack((np.zeros_like(heat_intl2),
                            heat_intl2)
                            ).reshape(heat_intl2.shape[0],-1)


    #FIGURE 1
    #figure configurations
    fig, ax = plt.subplots(figsize=[7,5])
    fig.patch.set_facecolor('white')

    #set the mask & labelling for the heatmap
    mask1 = np.vstack([np.arange(heat_intl1.shape[1])]* heat_intl1.shape[0]) % 3
    labelx = (np.repeat(expansionsteps,2)).tolist()[:mask1.shape[1]]
    labely = heatmapdata1.interaction_bp.unique().tolist()
    labely.sort()

    #plot settings
    plot_nr = sns.heatmap(heat_intl1,
                          mask = mask1,
                          ax = ax,
                          cmap=green,
                          linewidth=1,
                          cbar_kws={"shrink":0.4,'label': count_label},
                          annot_kws={"size": 8},
                          fmt="",
                          square=True,
                          vmin= 0.001,
                          #vmax= 100,
                          #vmax= heatmapdata1[count_type].max(),
                          xticklabels=labelx,
                          yticklabels=labely)

    plot_energy = sns.heatmap(heat_emin1,
                              mask = mask1-1,
                              ax = ax,
                              cmap=red,
                              linewidth=1, cbar_kws={"shrink":0.4, 'label': 'minE'},
                              annot_kws={"size": 8},
                              fmt="",
                              square=True,
                              vmin= heatmapdata1[energy_type].min(),
                              vmax= heatmapdata1[energy_type].max(),
                              xticklabels=labelx,
                              yticklabels=labely)

    plot_dif = sns.heatmap(heat_edif1,
                              mask = mask1-2,
                              ax = ax,
                              cmap=brown,
                              linewidth=1, cbar_kws={"shrink":0.4, 'label': 'minE'},
                              annot_kws={"size": 8},
                              fmt="",
                              square=True,
                              vmin= heatmapdata1[constraint_type].min(),
                              vmax= heatmapdata1[constraint_type].max(),
                              xticklabels=labelx,
                              yticklabels=labely)

    plt.title(title1, loc='center')
    plt.xlabel('constraint')
    plt.ylabel('interection length')

    plt.savefig(title1)


    #FIGURE 2: consecutive + nonconsecutive of the same run
    #figure configurations
    fig, ax = plt.subplots(figsize=[7,5])
    fig.patch.set_facecolor('white')

    shape1 = heat_intl12.shape
    shape2 = heat_intl2.shape
    a = shape1[0]-shape2[0]
    b = shape1[1]-shape2[1]

    def get_labels(heatmap):
        lx = list()
        c = 0

        for n in range(len(expansionsteps1)*2):
            if (n % 2 )== 0:
                lx.append(expansionsteps1[c])
                if expansionsteps1[c] == (heatmap["constraint"].max()):
                    lx.append(0)
                    lx.append(0)
                    break
                c += 1
            else:
                lx.append(0)

        ly = heatmap.interaction_bp.unique().tolist()
        ly.sort()

        return lx, ly


    if a < 0 or b < 0:
        heat_intl12 = np.pad(heat_intl12, ((0,a),(0,b)), 'constant', constant_values=(0))

        labelx, labely = get_labelx(heatmapdata2)

        mask2 = np.vstack([np.arange(heat_intl2.shape[1])]* heat_intl2.shape[0]) % 2

    elif a > 0 or b > 0:
        heat_intl2 = np.pad(heat_intl2, ((0,a),(0,b)), 'constant', constant_values=(0))

        labelx, labely = get_labelx(heatmapdata1)

        mask2 = np.vstack([np.arange(heat_intl12.shape[1])]* heat_intl12.shape[0]) % 2

    else:
        mask2 = np.vstack([np.arange(heat_intl12.shape[1])]* heat_intl12.shape[0]) % 2
        labelx = (np.repeat(expansionsteps,2)).tolist()[:mask2.shape[1]]
        labely = heatmapdata1.interaction_bp.unique().tolist()


    #plot settings
    plot_nr = sns.heatmap(heat_intl12,
                          mask = mask2,
                          ax = ax,
                          cmap=green,
                          linewidth=1,
                          cbar_kws={"shrink":0.4,'label': count_label},
                          annot_kws={"size": 8},
                          fmt="",
                          square=True,
                          vmin= 0.001,
                          #vmax= 100,
                          #vmax= heatmapdata1[count_type].max(),
                          xticklabels=labelx,
                          yticklabels=labely)

    plot_nr = sns.heatmap(heat_intl2,
                          mask = mask2-1,
                          ax = ax,
                          cmap=blue,
                          linewidth=1,
                          cbar_kws={"shrink":0.4,'label': count_label},
                          annot_kws={"size": 8},
                          fmt="",
                          square=True,
                          vmin= 0.001,
                          #vmax= 100,
                          #vmax= heatmapdata2[count_type].max(),
                          xticklabels=labelx,
                          yticklabels=labely)

    plt.title(title1, loc='center')
    plt.xlabel('constraint')
    plt.ylabel('interection length')

    plt.savefig(title2)



group_collect = list() #former fheatmap
constraint=list()

if clusters >= 1:
    # plot each cluster seperate
    ax1 = 0
    for step in stepsize1:
        new0_df= df_consecutive10[(df_consecutive10["stepsize"] == step)]

        title3 = (RNA+"_"+str(clusters)+consecutive+'_'+consecutive_run+"_"+str(step)+count_type)

        #figure setting
        fig, ax = plt.subplots(nrows=clusters,figsize=[12,40])
        fig.patch.set_facecolor('white')

        for cl in range(clusters):
            currentconstraint = expansionsteps[cl]

            new1_df= new0_df[(new0_df["cluster"] == cl)]

            for run in range(exsteps):
                new2_df = new1_df[(new1_df["constraint"] == expansionsteps[run])]

                group = new2_df.groupby(["interaction_countbp"]).agg({'interaction_countbp': ['count'],
                                                                    'energy_value': ['min', 'median', 'mean'],
                                                                    'energy_diff':['min','max','median','mean']}
                                                                 ).reset_index()
                group.columns =['interaction_bp', 'nr',
                                  'energy_min', 'energy_median','energy_mean',
                                  'dif_min','dif_max', 'dif_median','dif_mean']
                group['percent'] = (group.nr /group.nr.sum())*100
                group['constraint'] = expansionsteps[run]
                group.sort_values(by=['interaction_bp'])
                group=group.astype({'energy_min': float, 'energy_median': float,'energy_mean': float})
                group_collect.append(group)

            heatmapdata3 = pd.concat(group_collect)

            heat_intl3 = pd.pivot_table(heatmapdata3,
                                    index='interaction_bp',
                                    columns='constraint',
                                    values= count_type,
                                    aggfunc=np.sum,
                                    fill_value=0)

            heat_intl3 = np.dstack((heat_intl3,
                                    np.zeros_like(heat_intl3),
                                    np.zeros_like(heat_intl3))
                                ).reshape(heat_intl3.shape[0],-1)

            heat_emin3 = pd.pivot_table(heatmapdata3,
                                         index='interaction_bp',
                                         columns='constraint',
                                         values= energy_type,
                                         aggfunc=np.sum,
                                         fill_value=0)

            heat_emin3=np.dstack((np.zeros_like(heat_emin3),
                                    heat_emin3,
                                    np.zeros_like(heat_emin3))
                               ).reshape(heat_emin3.shape[0],-1)

            heat_edif3 = pd.pivot_table(heatmapdata3,
                                       index='interaction_bp',
                                       columns='constraint',
                                       values= constraint_type,
                                       aggfunc=np.sum,
                                       fill_value=0)

            heat_edif3=np.dstack((np.zeros_like(heat_edif3),
                                    np.zeros_like(heat_edif3),
                                    heat_edif3)
                                ).reshape(heat_edif3.shape[0],-1)


            #set the mask & labelling for the heatmap
            mask3 = np.vstack([np.arange(heat_intl3.shape[1])] * heat_intl3.shape[0]) % 3
            labelx = (np.repeat(expansionsteps,3)).tolist()[:mask3.shape[1]]
            labely = heatmapdata3.interaction_bp.unique().tolist()
            labely.sort()

            #FIGURE 1
            plot_nr = sns.heatmap(heat_intl3,
                                  mask = mask3,
                                  ax = ax[cl],
                                  cmap=green,
                                  linewidth=1,
                                  cbar_kws={"shrink":0.4,'label': count_label},
                                  annot_kws={"size": 8},
                                  fmt="",
                                  square=True,
                                  vmin= 0.001,
                                  vmax= heatmapdata3[count_type].max(),
                                  xticklabels=labelx,
                                  yticklabels=labely)

            plot_energy = sns.heatmap(heat_emin3,
                                      mask = mask3-1,
                                      ax = ax[cl],
                                      cmap=red,
                                      linewidth=1,
                                      cbar_kws={"shrink":0.4, 'label': 'minE'},
                                      annot_kws={"size": 8},
                                      fmt="",
                                      square=True,
                                      vmin= heatmapdata[energy_type].min(),
                                      vmax= heatmapdata[energy_type].max(),
                                      xticklabels=labelx,
                                      yticklabels=labely)

            plot_dif = sns.heatmap(heat_edif3,
                                  mask = mask3-2,
                                  ax = ax[cl],
                                  cmap=brown,
                                  linewidth=1,
                                  cbar_kws={"shrink":0.4, 'label': 'constraint'},
                                  annot_kws={"size": 8},
                                  fmt="",
                                  square=True,
                                  vmin= heatmapdata1[constraint_type].min(),
                                  vmax= heatmapdata1[constraint_type].max(),
                                  xticklabels=labelx,
                                  yticklabels=labely)


            title = title3
            plt.title(title, loc='center')
            plt.xlabel('constraint')
            plt.ylabel('interection length')
            group_collect.clear()


        plt.savefig(title3)
