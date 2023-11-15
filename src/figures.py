from functools import reduce

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

from flexitext import flexitext
from scipy.special import expit

import pyreadr

matplotlib.rcParams["font.family"] = "Times New Roman"
plasma_colormap = matplotlib.cm.get_cmap("plasma")

#COLORS = [plasma_colormap(x) for x in np.linspace(0.8, 0.15, num=4)]
#COLORS = [matplotlib.colors.to_hex(color) for color in COLORS]
COLORS = ["black", "blue", "blue"]

import sys

NOMPT = True
folder_name = "output"

try:
	isknown = sys.argv[1]
	sample_size = int(sys.argv[2])
	model_name = sys.argv[3]
except:
	isknown = "known"
	sample_size = 5000
	model_name = "m1"

if model_name == "m1":
	
	def f0(x, offset=-3):
		return np.sin(np.pi*x/3) + offset

	def f1(x, offset=-0.5):
		return np.pow(x, 3)/8 + offset

	def f2(x, offset=0):
		return -np.pow(x, 2)/4 + offset
	
if model_name == "m2":

	def f0(x, offset=-2.8):
		return np.sin(np.pi*x/3) + offset

	def f1(x, offset=0.5):
		return -np.pow(x,3)/8 + offset

	def f2(x, offset=-0.7):
		return np.pow(x, 2)/4 + offset

# This function is basically the code we wrote in the chunk above
def adjust_axis_layout(ax, title='', fn='f1', idx=[0,1], ntestings=3, sample_size=5000):
	# Remove major and minor tick marks on both axis
	ax.tick_params(axis="both", which="both", length=0)

	# Set major and minor ticks for the x axis.
	# These are used to draw the grid lines.
	# Only the major ticks have a tick label.
	ax.set_xticks([-2, -1, 0, 1, 2], minor=False)
	ax.set_xticklabels([-2, -1, 0, 1, 2], minor=False, size=11, color="0.3")
	ax.set_xticks([-2, -1, 0, 1, 2], minor=True)

	# Add grid lines for x axis
	ax.xaxis.grid(True, which="both", color="#cccccc", alpha=0.4, lw=0.5)

	# Set custom limit for x axis
	ax.set_xlim(-2.8, 2.8)

	if idx[1]==ntestings-1:
		twin = ax.twinx()
		twin.set_ylabel(f'$N={sample_size}$, {isknown} $S_e$ and $S_p$')
		twin.tick_params(left=False, right=False,
			labelleft=False, labelright = False) 
		for spine in ["top", "right", "bottom", "left"]:
			twin.spines[spine].set_visible(False)

	if fn == "f0":
		if model_name == "m1":
			# # Set custom limit for y axis
			ax.set_ylim(-5.6, -1)
			# Set major and minor ticks for the y axis.
			# The same logic than above.
			ax.set_yticks([-4.0, -3.0, -2.0], minor=False)
			ax.set_yticklabels([-4.0,  -3.0, -2.0], minor=False, size=11, color="0.3")
			ax.set_yticks([-4.0, -3.0, -2.0], minor=True)
		elif model_name == "m2":
			# # Set custom limit for y axis
			ax.set_ylim(-5.6, -1.0)
			# Set major and minor ticks for the y axis.
			# The same logic than above.
			ax.set_yticks([-4.0, -3.0, -2.0], minor=False)
			ax.set_yticklabels([-4.0,  -3.0, -2.0], minor=False, size=11, color="0.3")
			ax.set_yticks([-4.0, -3.0, -2.0], minor=True)

	elif fn == "f1":
		if model_name == 'm1':
			# Set custom limit for y axis
			#pass
			ax.set_ylim(-5, 3)
			# Set major and minor ticks for the y axis.
			# The same logic than above.
			ax.set_yticks([-4.0, -2.0, 0.0, 2.0], minor=False)
			ax.set_yticklabels([-4.0, -2.0, 0.0, 2.0], minor=False, size=11, color="0.3")
			ax.set_yticks([-4.0, -2.0, 0.0, 2.0], minor=True)
		elif model_name == 'm2':
			# Set custom limit for y axis
			#pass
			ax.set_ylim(-3, 4)
			# Set major and minor ticks for the y axis.
			# The same logic than above.
			ax.set_yticks([-2.0, 0.0, 2.0], minor=False)
			ax.set_yticklabels([-2.0, 0.0, 2.0], minor=False, size=11, color="0.3")
			ax.set_yticks([-2.0, 0.0, 2.0], minor=True)
	
	elif fn == "f2":
		if model_name == 'm1':
			# Set custom limit for y axis
			ax.set_ylim(-2.5, 0.4)
			# Set major and minor ticks for the y axis.
			# The same logic than above.
			ax.set_yticks([-2.0,-1.5,-1.0, -0.5, 0.0], minor=False)
			ax.set_yticklabels([-2.0,-1.5, -1.0, -0.5, 0.0], minor=False, size=11, color="0.3")
			ax.set_yticks([-2.0,-1.5,-1.0, -0.5, 0.0], minor=True)
		elif model_name == 'm2':
			# Set custom limit for y axis
			ax.set_ylim(-1.2, 1.5)
			# Set major and minor ticks for the y axis.
			# The same logic than above.
			ax.set_yticks([-0.5, 0.25, 1.0], minor=False)
			ax.set_yticklabels([-0.5, 0.25, 1.0], minor=False, size=11, color="0.3")
			ax.set_yticks([-0.5, 0.25, 1.0], minor=True)
	

	# Add grid lines for x axis
	ax.yaxis.grid(True, which="both", color="#cccccc", alpha=0.4, lw=0.5)

	# Remove all the spines
	for spine in ["top", "right", "bottom", "left"]:
		ax.spines[spine].set_visible(False)

	# if idx[1] == 0:
	# 	ax.spines["right"].set_visible(True)
	# if idx[1] == 1:
	# 	ax.spines["left"].set_visible(True)
	# if idx[1] == 2 and NOMPT:
	# 	ax.spines["right"].set_visible(True)
	# if idx[1] == 3 and NOMPT:
	# 	ax.spines["left"].set_visible(True)
	# if idx[1] == 3 and not NOMPT:
	# 	ax.spines["right"].set_visible(True)
	# if idx[1] == 4 and not NOMPT:
	# 	ax.spines["left"].set_visible(True)


	if idx[0] == 0:
		ax.set_title(title, weight=1000, size=12, loc="center")#, pad=0, y=1.000001)
	
	if idx[1] in range(1,ntestings,2):

		ax.set_facecolor(color="#f2f2f2")

	return ax

def get_estimates(df, label, testing, pool_size):

	df_reduce = df[(df["label"]==label)&(df["testing"]==testing)&(df["pool"]==pool_size)]
			
	knots = df_reduce["knots"]
	med = df_reduce["med"]
	lower = df_reduce["lower"]
	upper = df_reduce["upper"]
	truth = df_reduce["true"]

	return knots, truth, med, lower, upper





filename = f"{folder_name}/{isknown}/{sample_size}/{model_name}/df_pool.RData"

df = pyreadr.read_r(filename)[None]
# fig, ax = plt.subplots(figsize=(8, 6))

# knots = np.linspace(-2, 2, 200)
# discrete_fn_true = discrete_fn(knots)
# left = (knots < -1)
# right = (knots > 0.4)
# middle = ~(left | right)

# ax.plot(knots[left], discrete_fn_true[left], color=COLORS[0], alpha=1.0, lw=1.5)
# ax.plot(knots[middle], discrete_fn_true[middle], color=COLORS[0], alpha=1.0, lw=1.5)
# ax.plot(knots[right], discrete_fn_true[right], color=COLORS[0], alpha=1.0, lw=1.5)

# adjust_axis_layout(ax)



pool_sizes = [5, 10]
testings = ["IT", "DT", "AT", "MPT"]


if len(pool_sizes) == 1:
	testings = ["IT (c=1)", "DT (c=5)", "AT (c=5)", "MPT (c=5)"] if not NOMPT else ["IT (c=1)", "DT (c=5)", "AT (c=5)"]
else:
	testings = ["IT (c=1)", "DT (c=5)", "AT (c=5)", "MPT (c=5)", "DT (c=10)", "AT (c=10)", "MPT (c=10)"] if not NOMPT else ["IT (c=1)", "DT (c=5)", "AT (c=5)", "DT (c=10)", "AT (c=10)"]

labels = ["f0", "f1", "f2"]
lty = ["-","--", ":"]


# Initialize layout. Note we're using 1 row and 5 columns.
fig, axes = plt.subplots(len(labels), len(testings), 
	figsize=(16, 9), 
	sharey='row', 
	sharex=True)

# Set figure background color
fig.set_facecolor("white")

# Iterate over panels (programs)
# for j in range(6):
#	 # Select axis corresponding to the program
#	 ax = axes[j]
#	 # Create 100 replicates for each group
#	 for _ in range(100):
#		 probs = probabilities.compute(j)
#		 for prob, color in zip(probs, COLORS):
#			 ax.plot(x, prob, color=color, alpha=0.2, lw=1.2)
	
#	 # Note the title is unique for each panel/program
#	 adjust_axis_layout(ax, f"Program {j + 1}")



for i in range(len(labels)):
	
	label = labels[i]
	
	for j in range(len(testings)):
			
			testing = testings[j]
			
			if testing[:testing.index('(')-1] == 'IT':

				knots, truth, med, lower, upper = get_estimates(df, label, testing[:testing.index('(')-1], pool_size=10)
		
				ax = axes[i,j]

				ax.plot(knots, truth, color=COLORS[0], alpha=1.0, lw=1.5, linestyle=lty[0])
				ax.plot(knots, med, color=COLORS[0], alpha=1.0, lw=2, linestyle=lty[1])
				ax.fill_between(knots, lower, upper, color=COLORS[0], alpha=.1)

				ntestings = len(testings)
				adjust_axis_layout(ax, 
					title=testing, 
					fn=label, 
					idx=[i,j], 
					ntestings=ntestings,
					sample_size=sample_size)
				
				if i==2:
					ax.set_xlabel(f'$u$')
				if j==0:
					ax.set_ylabel(f'$\\beta_{i}$($u$)')

			else:
				pool = int(testing[(testing.index("=")+1):testing.index(")")])
				knots, truth, med, lower, upper = get_estimates(df, label, testing[:testing.index('(')-1], pool_size=pool)
		
				ax = axes[i,j]

				ax.plot(knots, truth, color=COLORS[0], alpha=1.0, lw=1.5, linestyle=lty[0])
				ax.plot(knots, med, color=COLORS[0], alpha=1.0, lw=2, linestyle=lty[1])
				ax.fill_between(knots, lower, upper, color=COLORS[0], alpha=.1)

				ntestings = len(testings)
				adjust_axis_layout(ax, 
					title=testing, 
					fn=label, 
					idx=[i,j], 
					ntestings=ntestings,
					sample_size=sample_size)
				
				if i==2:
					ax.set_xlabel(f'$u$')
				if j==0:
					ax.set_ylabel(f'$\\beta_{i}$($u$)')

			patch = mpatches.Patch(color=COLORS[0], alpha=0.2)

			line = [Line2D([0,1],[0,1],linestyle='-', color='black',lw=1.5),
				Line2D([0,1],[0,1],linestyle='--', color=COLORS[0], lw=2),
				patch]	

			if model_name=='m1':
				if i==0:
					location = "lower right"
				elif i==1:
					location = "lower right"
				else:
					location = "upper right"
			if model_name=='m2':
				if i==0:
					location = "upper right"
				elif i==1:
					location = "upper right"
				else:
					location = "lower right"
			if i==0:
				ax.legend(handles=line,
					labels=['Truth', 'Median', r"95\% CI"], 
					ncol=3 if NOMPT else 1, 
					fontsize="8",
					frameon=False,
					# borderpad=1.3,
					# columnspacing=3,
					loc="upper center")

fig.tight_layout()

#fig.subplots_adjust(top=0.85)

if NOMPT:
	plt.savefig(f"{isknown}_{model_name}_wo_MPT.png", format="png", dpi=300)
else:
	plt.savefig(f"{isknown}_{model_name}_w_MPT.png", format="png", dpi=300)

plt.show()





