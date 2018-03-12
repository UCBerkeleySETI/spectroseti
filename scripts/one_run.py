import matplotlib
matplotlib.use('TkAgg')

import spectroseti.runner as runner


observations =[['baa',201],['baa',202],['baa',203],['baa',204]]#,['awx',222],['awx',223],['awx',224]]


LS = runner.LaserSearch()

LS.search_multiple(observations,output_pngs=1)