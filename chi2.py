import numpy as np
import matplotlib.pyplot as plt
import sys
#from matplotlib.ticker import FormatStrFormatter
#from scipy.stats import norm
#import matplotlib.mlab as mlab
#import csv
import pandas as pd

infile = sys.argv[1]
infile2 = sys.argv[2]
#infile3 = sys.argv[3]

data = pd.read_csv(infile, sep=",")
data2= pd.read_csv(infile2, sep=",")
#data3= pd.read_csv(infile3, sep=",")

temp_list=[]
for i in range(0,16):
       temp_list.append("data_channel{0}".format(i))

for x in range(0, 8):
    df=pd.DataFrame(data)
    df2=pd.DataFrame(data2)
    #df3=pd.DataFrame(data3)

    df = df[df.channel == x]
    df2 = df2[df2.channel == x]
    #df2 = df2[df2.ch == x]

    #df3 = df3[df3.channel == x]

    plt.figure();
    #plt.ylim(0.6,2.4)
    plt.plot(df["freqs"], abs(df["p1"]), label=temp_list[x]) 
    #    pfit = np.polyfit(df["freqs"], abs(df["p1"])/abs(df2["p1"]), 0)
    # plt.plot(df["freqs"], abs(df2["p1"]), label=temp_list[x])
    #  trend_line_model = np.poly1d(pfit)
    #  plt.plot(df["freqs"], trend_line_model(df["freqs"]), "m--")
    # print(np.mean(abs(df["p1"])/abs(df2["p1"])))
    plt.plot(df2["freqs"], abs(df2["p1"]),"--", label='sigmas_from_AraSim')
   # plt.plot(df2["Freq"], abs(df2["sigma"]),"--", label='sigmas_from_AraSim')
   #plt.plot(df2["freqs"], abs(df2["p1"])/abs(df["p1"]),"", label='Quotient AraSim/Data')

    #plt.plot(df3["freqs"], df3["p1"],label='AraSim_jorge')
  # print(np.mean(abs(df2["p1"])/abs(df["p1"])))
    #print(df["p1"]/df2["p1"])
    plt.ylabel("$\sigma$ of Rayleigh dist.")
    plt.xlabel("Frequency [MHz]")
    #plt.ylim(0.,0.5)
    plt.legend(loc='best')
plt.title("Freq vs $\sigma$")
plt.savefig("sigmavsfreq_ch0.png",bbox_inches="tight")
plt.show()

