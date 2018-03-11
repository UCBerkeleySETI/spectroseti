import matplotlib
matplotlib.use('TkAgg')

import spectroseti.apf as apf
import spectroseti.utilities as util
import random
import numpy as np

import matplotlib.pyplot as plt
import pandas as pd
import time
import Tkinter
import tkMessageBox



y = 15
x = 10


obs = [['awx',219],['amp',193],['ayo',224]]
# Do this for three observations:

crs = pd.DataFrame()
fig = plt.figure()

cosmic = np.expand_dims(np.zeros((y*2,x*2)),2)
no_cosmic = np.expand_dims(np.zeros((y*2,x*2)),2)

for ob in obs:
    raw = apf.APFRawObs('awx', 219)
    print(' New Star ')
    for i in range(100):
        # Pick a random location
        subset = raw.retrieve_subset(random.randint(0,78), random.randint(30,4577),yradius=y,xradius=x)
        # plot a nxn image

        #render and await input 'y', 'n' (y has CR, n has no CR)
        #root=Tkinter.Tk()
        plt.imshow(subset)
        fig.canvas.flush_events()
        cb = plt.colorbar()

        plt.show(block=False)

        fig.canvas.flush_events()
        # Option messagebox
        response = raw_input("<Hit Enter To Close>")
        #response = tkMessageBox.askyesno("CR Rejector", "Is there a cosmic ray?")
        #print(response)
        #root.destroy()
        plt.cla()
        cb.remove()

        # There IS a cosmic
        if response == '1':
            print('cosmic')
            cosmic = np.append(cosmic, np.expand_dims(subset, 2), axis=2)
            print cosmic.shape

        # there is NO cosmic
        elif response == '2':
            print('no cosmic')
            print no_cosmic.shape, subset.shape
            no_cosmic = np.append(no_cosmic, np.expand_dims(subset, 2), axis=2)
            print no_cosmic.shape


np.save('cosmic_examples.npy', cosmic)
np.save('no_cosmic_examples.npy', no_cosmic)

# accumulate 12x6 ? postage stamps in "cosmic", "noncosmic" dataframes


# save these as .npy
