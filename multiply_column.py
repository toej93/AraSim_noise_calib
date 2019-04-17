# -*- coding: utf-8 -*-
import numpy as np
import sys
#import matplotlib.pyplot as plt
from pylab import setp
from matplotlib.pyplot import rcParams
import csv
import glob
import pandas as pd

data = pd.read_csv("Rayleigh_A2.csv", sep=",")
df=pd.DataFrame(data)
df.loc[df.channel == 15,'p1']*=0.001
#print(df.head(5))
#df.loc[:,'p1'] *=0.001 
print(df.loc[df.channel == 15,'p1'].head(5))
df.to_csv('Rayleigh_A2_masked.csv', sep=',',index=False)

    
