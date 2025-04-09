import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set font to Times New Roman
rcParams['text.usetex'] = True  # Enable LaTeX rendering
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']

# Load data
huizong_A= np.loadtxt('trajectory_A_cells_0.25')
huizong_B = np.loadtxt('trajectory_B_cells_0.25')
#huizong = np.loadtxt('huizong_0.25')
huizong1 = np.loadtxt('huizong_A_0.25')
huizong2 = np.loadtxt('huizong_B_0.25')

# Create time points
time_points = np.arange(0, 39.1, 0.1)

# Create high-quality figure
fig = plt.figure(figsize=(10,5),dpi=300)
ax1= fig.add_subplot(1,2,1)

# Plot all trajectories with low alpha for transparency
for j in range(huizong_A.shape[1])
    ax1.plot(time_points,  np.log(huizong_A[, j]), 'c-', linewidth=0.5, alpha=0.1)

for j in range(huizong_B.shape[1])
    ax1.plot(time_points, np.log( huizong_B[, j]), 'm-', linewidth=0.5, alpha=0.1)

# Plot mean lines with thicker linewidth
mean_1 = np.mean(huizong_A,axis=1)
mean_2 = np.mean(huizong_B, axis=1)

ax1.plot(time_points, np.log( mean_1), 'c-', linewidth=1, label=rAverage $X_t$)
ax1.plot(time_points, np.log(mean_2), 'm-', linewidth=1, label=rAverage $Y_t$)
# 将对数刻度转换为原始值显示在 y 轴上

yticks = ax1.get_yticks()  # 获取当前的 y 轴刻度
ax1.set_yticklabels([f'{np.exp(tick).2e}' for tick in yticks])  # 将对数值转换为原始值
# Add grid
ax1.grid(True, linestyle='--', alpha=0.5)

# Set labels with LaTeX formatting
ax1.set_xlabel('time $t$', fontsize=16)
ax1.set_ylabel(r'cell count',fontsize=16)

# Set tick parameters
#ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.set_title((c),fontsize=16)
# Add legend with larger font
ax1.legend(fontsize=16,loc = upper left)

# Adjust layout
plt.tight_layout()

ax2= fig.add_subplot(1,2,2)

# Plot all trajectories with low alpha for transparency
for j in range(huizong1.shape[1])
    ax2.plot(time_points, huizong1[, j], 'c-', linewidth=0.5, alpha=0.1)

for j in range(huizong2.shape[1])
    ax2.plot(time_points, huizong2[, j], 'm-', linewidth=0.5, alpha=0.1)

# Plot mean lines with thicker linewidth
mean_1 = np.mean(huizong1, axis=1)
mean_2 = np.mean(huizong2, axis=1)

ax2.plot(time_points, mean_1, 'c-', linewidth=1, label=r Average $mu_t$)
ax2.plot(time_points, mean_2, 'm-', linewidth=1, label=rAverage $gamma_t$)

# Add grid
ax2.grid(True, linestyle='--', alpha=0.5)

# Set labels with LaTeX formatting
ax2.set_xlabel('time $t$', fontsize=16)
ax2.set_ylabel(r'cell proportion', fontsize=16)
ax2.set_title((d),fontsize=16)
# Set tick parameters
#ax2.tick_params(axis='both', which='major', labelsize=14)

# Add legend with larger font
ax2.legend(fontsize=16,loc = upper left)

# Optional Add text annotations if needed
# ax.text(25, max(mean_1)0.8, r$alpha = 0.35$, fontsize=22, ha='center')

# Adjust layout
plt.tight_layout()
plt.savefig( 'result2_traj.pdf',dpi=300)#第一个是指存储路径，第二个是图片名字
plt.show()