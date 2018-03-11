import matplotlib
matplotlib.use('TkAgg')

import spectroseti.runner as runner


observations =[['awx',221]]#,['awx',222],['awx',223],['awx',224]]


LS = runner.LaserSearch()

LS.search_multiple(observations,output_pngs=1)