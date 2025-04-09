import numpy as np
import math
import matplotlib.pyplot as plt

# Load data
huizong1 = np.loadtxt('huizong_A_0.25')
huizong2 = np.loadtxt('huizong_B_0.25')
huizong = np.loadtxt('huizong_0.25')

# Second plot for huizong with mean and standard deviation
time_points = np.arange(0, 39.1, 0.1)

# First plot for huizong1 and huizong2
average_1 = []
average_2 = []

for h in range(len(huizong1)):
    average1 = np.mean(huizong1[h])
    average_1.append(average1)

for i in range(len(huizong2)):
    average2 = np.mean(huizong2[i])
    average_2.append(average2)
Nt_breed1 = []
std_values = []
for h in range(len(huizong)):
    mean1 = np.mean(huizong[h])
    std1 = np.std(huizong[h])
    Nt_breed1.append(mean1)
    std_values.append(std1)

# Convert to numpy arrays for easier manipulation
Nt_breed1 = np.array(Nt_breed1)
std_values = np.array(std_values)

n = 1
ave_new = []
time_new = []
tt = []
for i in range(0, len(average_1), n):
    ave_1 = average_1[i:i + n]
    ave_new.append(ave_1)
for j in np.arange(0, 39.1, 0.1):
    tt.append(j)
for k in range(0, len(tt), n):
    time_1 = tt[k:k + n]
    time_new.append(time_1)
var_1=[]
var_2=[]
Nt_var1 = []
lambda1 = 0.4
lambda2 = 0.2
alpha = 0.35
beta = 0.25
for h in range(0,len(huizong1)):
    average1= np.mean(huizong1[h])
    var1= np.var(huizong1[h])
    Nt_var_1 = ((2*lambda1*alpha-2*lambda2*beta+lambda2-lambda1)*average1-lambda2+2*lambda2*beta)/(2*(lambda2-lambda1)*var1)
    Nt_var1.append(Nt_var_1)

n = 1
ave_new = []
time_new = []
tt = []
for i in range(0, len(average_1), n):
    ave_1 = average_1[i:i + n]
    ave_new.append(ave_1)
for j in np.arange(0, 39.1, 0.1):
    tt.append(j)
for k in range(0, len(tt), n):
    time_1 = tt[k:k + n]
    time_new.append(time_1)

m = len(ave_new)-1
Nt_mean1 = np.zeros(m)
N0 = 100
delta_t = 0.1*n
sum_mt = []
Nt_mean1[0] = N0
for k in range(1, m):
    mt_1 = []
    mt_n = 0
    for i in range(0, k):
        mt_1.append((ave_new[i][0] + ave_new[i + 1][0]) / 2)
    for g in range(0, len(mt_1)):
        mt_n = mt_n + mt_1[g]
    sum_mt.append(mt_n)
    Nt_mean1[k] = N0*np.exp((lambda1-lambda2)*delta_t*sum_mt[-1]+lambda2*time_new[k][0])

# 设置字体为罗马字体（Times New Roman）
rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']

# 创建图形和子图
fig = plt.figure(figsize=(10,5), dpi=300)

# 画第1个图：折线图
ax1 = fig.add_subplot(121)
ax1.plot(np.arange(0, 39, 0.1), np.log(Nt_mean1), c = "c", label = "FOM")
ax1.plot(np.arange(0, 39.1, 0.1), np.log(Nt_breed1), c = "k", linewidth=1, linestyle="--", label="Simulation")


ax1.fill_between(time_points, 
                 np.log(Nt_breed1 - std_values), 
                 np.log(Nt_breed1 + std_values),
                 color="gray", alpha=0.3, label="±SD")

#ax1.text(25, 6, "$MSE_{F}$ :  2e-05", fontsize=14, ha='center')
ax1.text(20, 7, "$MRE_{F}$ :  0.01203", fontsize=14, ha='left')
ax1.text(20, 6, "$Slope_{F}$ :  0.261", fontsize=14, ha='left')
ax1.text(18, 5, "$Slope_{Sim}$ :  0.261", fontsize=14, ha='left')
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.set_xlabel('time t', fontsize=14, color='k')
ax1.set_ylabel('$N_t$', fontsize=14, color='k')
ax1.set_title("(a)", fontsize=14)

# 将对数刻度转换为原始值显示在 y 轴上
yticks = ax1.get_yticks()  # 获取当前的 y 轴刻度
ax1.set_yticklabels([f'{np.exp(tick):.2e}' for tick in yticks])  # 将对数值转换为原始值

plt.legend(fontsize=12)  

# 画第2个图：散点图
ax2 = fig.add_subplot(122)
ax2.plot(np.arange(0, 39.1, 0.1), np.log(Nt_breed1), c = "k", linewidth=1, linestyle="--", label="Simulation")
#ax2.plot(np.arange(0, 39.1, 0.1), np.log(Nt_breed1), c = "k", linewidth=1, linestyle="--", label="Simulation")

ax2.fill_between(time_points, 
                 np.log(Nt_breed1 - std_values), 
                 np.log(Nt_breed1 + std_values),
                 color="gray", alpha=0.3, label="±SD")
ax2.plot(np.arange(0, 39.1, 0.1), np.log(Nt_var1), c = "red", label="SOM")
#ax2.text(25, 6, "$MSE_{S}$ :  0.02917", fontsize=14, ha='center')
ax2.text(20, 7, "$MRE_{S}$ :  0.13062", fontsize=14, ha='left')
ax2.text(20, 6, "$Slope_{S}$ :  0.253", fontsize=14, ha='left')
ax2.text(18, 5, "$Slope_{Sim}$ :  0.261", fontsize=14, ha='left')
ax2.grid(True, linestyle='--', alpha=0.5)
ax2.set_xlabel('time t', fontsize=14, color='k')
ax2.set_title("(b)", fontsize=14)

# 将对数刻度转换为原始值显示在 y 轴上
yticks2 = ax2.get_yticks()  # 获取当前的 y 轴刻度
ax2.set_yticklabels([f'{np.exp(tick):.2e}' for tick in yticks2])  # 将对数值转换为原始值

plt.legend(fontsize=12)  
#plt.savefig( 'model1_result2_SD.pdf',dpi=300)#第一个是指存储路径，第二个是图片名字
# 显示图像
plt.show()