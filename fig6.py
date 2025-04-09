import numpy as np
import math
import matplotlib.pyplot as plt
huizong1=np.loadtxt('huizong_A_0.25')
huizong2=np.loadtxt('huizong_B_0.25')
huizong=np.loadtxt('huizong_0.25')

average_1=[]
average_2=[]
for h in range(len(huizong1)):
    average1= np.mean(huizong1[h])
    average_1.append(average1)

for i in range(len(huizong2)):
    average2= np.mean(huizong2[i])
    average_2.append(average2)
Nt_breed1 = []
std_values = []
time_points = np.arange(0, 39.1, 0.1)
for h in range(len(huizong)):
    mean1 = np.mean(huizong[h])
    std1 = np.std(huizong[h])
    Nt_breed1.append(mean1)
    std_values.append(std1)

# Convert to numpy arrays for easier manipulation
Nt_breed1 = np.array(Nt_breed1)
std_values = np.array(std_values)

# 创建新列表来保存添加了噪声的数据
huizong1_noisy = []
huizong2_noisy = []

# 遍历并为每个子列表中的数据添加噪声
for i in range(len(huizong1)):
    # 创建新的子列表来保存添加了噪声的数据
    noisy_sublist = []
    for j in range(len(huizong1[i])):
        noise = np.random.normal(0, 0.01)
        noisy_value = huizong1[i][j] + noise
        noisy_sublist.append(noisy_value)
    huizong1_noisy.append(noisy_sublist)

for i in range(len(huizong2)):
    # 创建新的子列表来保存添加了噪声的数据
    noisy_sublist = []
    for j in range(len(huizong2[i])):
        noise = np.random.normal(0, 0.01)
        noisy_value = huizong2[i][j] + noise
        noisy_sublist.append(noisy_value)
    huizong2_noisy.append(noisy_sublist)
#时间间隔增大为原来的倍数

n = 1
ave_new = []
time_new = []
tt = []
for i in range(0, len(average_n1), n):
    ave_1 = average_n1[i:i + n]
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
for h in range(0,len(huizong1_noisy),n):
    average1= np.mean(huizong1_noisy[h])
    var1= np.var(huizong1_noisy[h])
    Nt_var_1 = ((2*lambda1*alpha-2*lambda2*beta+lambda2-lambda1)*average1-lambda2+2*lambda2*beta)/(2*(lambda2-lambda1)*var1)
    Nt_var1.append(Nt_var_1)

n = 1
ave_new = []
time_new = []
tt = []

# 对 average_1 进行处理
for i in range(0, len(average_n1), n):
    ave_1 = average_1[i:i + n]
    ave_new.append(ave_1)

# 生成时间序列
for j in np.arange(0, 39.1, 0.1):
    tt.append(j)

for k in range(0, len(tt), n):
    time_1 = tt[k:k + n]
    time_new.append(time_1)

# 计算 Nt_mean1
m = len(ave_new) - 1
Nt_mean1 = np.zeros(m)
N0 = 100
delta_t = 0.1 * n  # 调整时间间隔

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
    Nt_mean1[k] = N0 * np.exp((lambda1 - lambda2) * delta_t * sum_mt[-1] + lambda2 * time_new[k][0])

#现在要在一张图上对数拟合曲线，然后再图中标出MSE，再标出平均相对误差
#尝试用MSE来评估两条曲线的拟合程度
def calculate_mse(data_points, fitted_curve):
    # 将列表转换为 NumPy 数组
    data_points = np.array(data_points)
    fitted_curve = np.array(fitted_curve)
    
    # 计算差值
    differences = data_points - fitted_curve
    # 计算差值的平方
    squared_differences = differences ** 2
    # 计算均方误差
    mse = np.mean(squared_differences)
    return mse

#用平均相对误差来衡量曲线的拟合程度
# 计算相对误差
def calculate_relative_error_per_point(data_points2, fitted_curve2):
    # 将列表转换为 NumPy 数组
    data_points2 = np.array(data_points2)
    fitted_curve2 = np.array(fitted_curve2)
    
    # 计算绝对误差
    absolute_error2 = np.abs(data_points2 - fitted_curve2)
    
    # 计算相对误差
    relative_error= absolute_error2 / np.abs(data_points2)
    
    return relative_error
# 示例数据点和拟合曲线
data_points =np.log(Nt_breed1[:-1])
fitted_curve_v = np.log(Nt_var1[:-1])
fitted_curve_m =np.log(Nt_mean1) 

data_points2 =Nt_breed1[:-1]
fitted_curve_v2 = Nt_var1[:-1]
fitted_curve_m2 =Nt_mean1

# 计算均方误差
mse_v = calculate_mse(data_points, fitted_curve_v)
mse_m = calculate_mse(data_points, fitted_curve_m)
# 计算相对误差
re_v = calculate_relative_error_per_point(data_points2, fitted_curve_v2)
re_m = calculate_relative_error_per_point(data_points2, fitted_curve_m2[:len(data_points)])  # 确保长度一致

# 计算平均相对误差
mre_v = np.mean(re_v)
mre_m = np.mean(re_m)

rounded_mse_v = round(mse_v,5)
rounded_mse_m = round(mse_m,5)
rounded_mre_v = round(mre_v,5)
rounded_mre_m = round(mre_m,5)

from matplotlib import rcParams

# 设置字体为罗马字体（Times New Roman）
rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']

fig = plt.figure(figsize=(6,5),dpi=300)
ax = fig.add_subplot()  
ax.plot(np.arange(0,39.1,0.1),np.log(Nt_breed1), c = "k",linewidth=3,linestyle="--",label="Simulation") 
ax.plot(np.arange(0,39.1,0.1),np.log(Nt_var1), c = "red",label="SOM")
ax.plot(np.arange(0,39,0.1),np.log(Nt_mean1),c = "c",label = "FOM")
ax.text(25,5,"$MRE_{S}$ : 0.68826 ",fontsize=22, ha='center')
ax.text(25,6,"$MRE_{F}$ : 0.00443",fontsize=22, ha='center')
#plt.gca().set(xlabel='time t', ylabel= 'log$N_t$')
#plt.title("noise ~ $N(0,0.01)$",fontsize=14)
ax.legend(fontsize=22)
ax.set_xlabel('time t',fontsize=24)
ax.set_ylabel('$N_t$',fontsize=24)
ax.grid(True, linestyle='--', alpha=0.5)

yticks = ax.get_yticks()  # 获取当前的 y 轴刻度
ax.set_yticklabels([f'{np.exp(tick):.2e}' for tick in yticks])  # 将对数值转换为原始值
plt.savefig( 'noise(MRE)_1.pdf',dpi=300,bbox_inches='tight')
plt.show()