import matplotlib.pyplot as plt
import numpy as np
import re
import matplotlib as mpl

# 设置全局字体大小
mpl.rcParams.update({'font.size': 14})

round = [1,2,3]
num_x = [1,2,3,4,5,6,7,8,9,10]
verifier_time_deg_1 = []
verifier_time_deg_2 = []
verifier_time_deg_4 = []

def extract_number(line):
    # 使用正则表达式提取数字
    pattern = r"Sumcheck on F\+G Verifier time: ([0-9.]+)"
    match = re.search(pattern, line)

    if match:
        time_data = float(match.group(1))
        return time_data
    else:
        return None


with open('experiments.txt', 'r') as file:
    content = file.readlines()

current_experiment = None
#获取experment 1-1的prover time
for line in content:  
    if line.startswith('-------------experiment 1-1'):
         current_experiment = "experiment 1-1"
    if line.startswith('Sumcheck on F+G Verifier time') and current_experiment == "experiment 1-1":             
            time_data = extract_number(line)
            verifier_time_deg_1.append(time_data)         
    elif line.startswith('-------------experiment 1-2'):
        break

#获取experment 1-2的prover time
for line in content:  
    if line.startswith('-------------experiment 1-2'):        
        current_experiment = "experiment 1-2"
    elif line.startswith('Sumcheck on F+G Verifier time') and current_experiment == "experiment 1-2":             
            time_data = extract_number(line)
            verifier_time_deg_2.append(time_data)         
    elif line.startswith('-------------experiment 1-3'):
        break
#获取experment 1-3的prover time
for line in content:  
    if line.startswith('-------------experiment 1-3'):
         current_experiment = "experiment 1-3"
    elif line.startswith('Sumcheck on F+G Verifier time') and current_experiment == "experiment 1-3":             
            time_data = extract_number(line)
            verifier_time_deg_4.append(time_data)         
    elif line.startswith('-------------experiment 2'):
        break
                



plt.figure(figsize=(4, 3))

 # 画折线图
plt.plot(num_x, verifier_time_deg_1,  marker='o', linestyle='-')
plt.plot(num_x, verifier_time_deg_2,  marker='o', linestyle='-')
plt.plot(num_x, verifier_time_deg_4,  marker='o', linestyle='-')

selected_ticks = [2, 4, 6, 8, 10]
plt.xticks(selected_ticks, ['$2^2$', '$2^4$', '$2^6$', '$2^8$', '$2^{10}$'])
#plt.xticks(num_x, ['$2^1$', '$2^2$', '$2^3$', '$2^4$', '$2^5$', '$2^6$', '$2^7$', '$2^8$', '$2^9$','$2^{10}$'])
plt.xlabel('k')
plt.ylabel('time (ms)')
plt.legend()
plt.legend(frameon=False)
plt.tight_layout()
# 显示图形
plt.grid(True)
plt.show()
