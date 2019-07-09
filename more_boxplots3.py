from matplotlib.pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
import json

import matplotlib.pyplot as plt


def _read_file(name, iteration):
    arr_of_auc = []
    with open(name + "/" + iteration + "_analysis/auc.json", "r", encoding="utf-8") as input_stream:
            for new_line in input_stream:
                line = json.loads(new_line)
                arr_of_auc.append(line["AUC"])
    return arr_of_auc                

                
# function for setting the colors of the box plots pairs
def setBoxColors(bp, ok: bool):
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['fliers'][0], markeredgecolor='blue',markersize=0.15)
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    setp(bp['fliers'][1], markeredgecolor='red', markersize=0.15)
    setp(bp['medians'][1], color='red')
    
    setp(bp['boxes'][2], color='brown')
    setp(bp['caps'][4], color='brown')
    setp(bp['caps'][5], color='brown')
    setp(bp['whiskers'][4], color='brown')
    setp(bp['whiskers'][5], color='brown')
    setp(bp['fliers'][2], markeredgecolor='brown', markersize=0.15)
    setp(bp['medians'][2], color='brown')
    
    setp(bp['boxes'][3], color='green')
    setp(bp['caps'][6], color='green')
    setp(bp['caps'][7], color='green')
    setp(bp['whiskers'][6], color='green')
    setp(bp['whiskers'][7], color='green')
    setp(bp['fliers'][3], markeredgecolor='green', markersize=0.15)
    setp(bp['medians'][3], color='green')

# Some fake data to plot
A= [[2],  [2], [2], [2]]
B = [[5, 7, 2, 2, 5], [7, 2, 5], [7,3,4], [2,3]]
C = [[3,2,5,7], [6, 7, 3], [5,8,8], [4,5]]
D= [[1, 2, 5],  [7, 2], [8,9], [1,4]]
E = [[5, 7, 2, 2, 5], [7, 2, 5], [7,3,4], [2,3]]
F = [[3,2,5,7], [6, 7, 3], [5,8,8], [4,5]]
G= [[1, 2, 5],  [7, 2], [8,9], [1,4]]
H = [[5, 7, 2, 2, 5], [7, 2, 5], [7,3,4], [2,3]]
I = [[3,2,5,7], [6, 7, 3], [5,8,8]]

fig = figure()
fig.set_size_inches(6, 3)
ax = axes()

A1 = [[0.9267435897436331]]
bp = boxplot(A1, positions=[0.1], widths= 0.19)

# first boxplot pair
A = [[0.9178333333333628], [0.9178333333333628], [0.9236282051282391], [0.9236282051282391]] 
bp = boxplot(A, positions=[0.7,0.9,1.1,1.3], widths= 0.19)
setBoxColors(bp, True)

# second boxplot pair
iter21 = _read_file("cutoff_active_group_model0_e4", "pairs")
iter22 = _read_file("cutoff_and_group_model0", "pairs")
iter23 = _read_file("cutoff_active_group_model", "pairs")
iter24 = _read_file("cutoff_and_group_model_100", "pairs")
B = [iter21, iter22, iter23, iter24]

bp = boxplot(B, positions=[1.7,1.9,2.1,2.3], widths= 0.19)
setBoxColors(bp, True)

# thrid boxplot pair
iter21 = _read_file("cutoff_active_group_model0_e4", "triples")
iter22 = _read_file("cutoff_and_group_model0", "triples")
iter23 = _read_file("cutoff_active_group_model", "triples")
iter24 = _read_file("cutoff_and_group_model_100", "triples")
C = [iter21, iter22, iter23, iter24]

bp = boxplot(C, positions=[2.7,2.9,3.1,3.3], widths= 0.19)
setBoxColors(bp, True)

iter21 = _read_file("cutoff_active_group_model0_e4", "quatres")
iter22 = _read_file("cutoff_and_group_model0", "quatres")
iter23 = _read_file("cutoff_active_group_model", "quatres")
iter24 = _read_file("cutoff_and_group_model_100", "quatres")
D = [iter21, iter22, iter23, iter24]

bp = boxplot(D, positions=[3.7,3.9,4.1,4.3], widths= 0.19)
setBoxColors(bp, True)

iter21 = _read_file("cutoff_active_group_model0_e4", "fifth")
iter22 = _read_file("cutoff_and_group_model0", "fifth")
iter23 = _read_file("cutoff_active_group_model", "fifth")
iter24 = _read_file("cutoff_and_group_model_100", "fifth")
E = [iter21, iter22, iter23, iter24]

bp = boxplot(E, positions=[4.7,4.9,5.1,5.3], widths= 0.19)
setBoxColors(bp, True)

iter21 = _read_file("cutoff_active_group_model0_e4", "sixth")
iter22 = _read_file("cutoff_and_group_model0", "sixth")
iter23 = _read_file("cutoff_active_group_model", "sixth")
iter24 = _read_file("cutoff_and_group_model_100", "sixth")
F = [iter21, iter22, iter23, iter24]

bp = boxplot(F, positions=[5.7,5.9,6.1,6.3], widths= 0.19)
setBoxColors(bp, True)

iter21 = _read_file("cutoff_active_group_model0_e4", "septh")
iter22 = _read_file("cutoff_and_group_model0", "septh")
iter23 = _read_file("cutoff_active_group_model", "septh")
iter24 = _read_file("cutoff_and_group_model_100", "septh")
G = [iter21, iter22, iter23, iter24]

bp = boxplot(G, positions=[6.7,6.9,7.1,7.3], widths= 0.19)
setBoxColors(bp, True)

iter21 = _read_file("cutoff_active_group_model0_e4", "octagon")
iter22 = _read_file("cutoff_and_group_model0", "octagon")
iter23 = _read_file("cutoff_active_group_model", "octagon")
iter24 = _read_file("cutoff_and_group_model_100", "octagon")
H = [iter21,  iter22, iter23, iter24]

bp = boxplot(H, positions=[7.7,7.9,8.1,8.3], widths= 0.19)
setBoxColors(bp, True)

iter21 = _read_file("cutoff_active_group_model0_e4", "ninth")
iter22 = _read_file("cutoff_and_group_model0", "ninth")
iter23 = _read_file("cutoff_active_group_model", "ninth")
iter24 = _read_file("cutoff_and_group_model_100", "ninth")
I = [iter21, iter22, iter23, iter24]

bp = boxplot(I, positions=[8.7,8.9,9.1, 9.3], widths = 0.19)
setBoxColors(bp, True)

# set axes limits and labels

ax.set_xticklabels(['', '1.', '2.', '3.', '4.', '5.', '6.', '7.','8.','9.'])
ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
plt.axvline(x=0.5, color="black", linewidth=0.5)
plt.axvline(x=1.5, color="black", linewidth=0.5)
plt.axvline(x=2.5, color="black", linewidth=0.5)
plt.axvline(x=3.5, color="black", linewidth=0.5)
plt.axvline(x=4.5, color="black", linewidth=0.5)
plt.axvline(x=5.5, color="black", linewidth=0.5)
plt.axvline(x=6.5, color="black", linewidth=0.5)
plt.axvline(x=7.5, color="black", linewidth=0.5)
plt.axvline(x=8.5, color="black", linewidth=0.5)

# draw temporary red and blue lines and use them to create a legend
#hB, = plot([1,1],'b-')
#hR, = plot([1,1],'r-')
#hB.set_visible(False)
#R.set_visible(False)
ax.legend([bp["boxes"][0],bp["boxes"][1],bp["boxes"][2],bp["boxes"][3]], ["O_0", "OP_0", "O_40", "OP_40"], loc="upper left", prop={'size': 8})
plt.ylabel("AUC")
plt.xlabel("iterace")
plt.tight_layout()
savefig('boxcompare3.png', dpi=1000)
show()