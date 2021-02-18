import matplotlib.pyplot as plt

def stats_plot(entangled,psc,protid,fileformat='jpeg'):
 
    fig,(ax1,ax2)= plt.subplots(1,2)
    ax1.plot(entangled)
    ax2.pie(psc,autopct = '%1.1f%%')
    ax1.set_xlabel('Distance from diagonal')
    ax1.set_ylabel('Fraction entangled')
    fig.suptitle(protid)
    ax2.legend(['Paralell','Series','Cross'],bbox_to_anchor=(.5, -0.5, 0.5, 0.5) )

