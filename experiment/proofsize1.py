import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

mpl.rcParams.update({'font.size': 14})
# 数据
x_values = np.arange(1, 10)
d_values = [1,2,4] # degree of h
log_k_values = 2  # k instances
log_m_values = np.arange(1, 10)  # m-length witness
m_values = 2 ** log_m_values
F_size = 256/8/1024 # KB
G_size = 512/8/1024 # KB

# 优化证明大小
optimized_proof_size_F_d1 = ((d_values[0] + 1) * log_k_values + 7)
optimized_proof_size_G_d1 = ((d_values[0] + 1) * log_m_values + 2)

optimized_proof_size_F_d2 = ((d_values[1] + 1) * log_k_values + 7)
optimized_proof_size_G_d2 = ((d_values[1] + 1) * log_m_values + 2)

optimized_proof_size_F_d4 = ((d_values[2] + 1) * log_k_values + 7)
optimized_proof_size_G_d4 = ((d_values[2] + 1) * log_m_values + 2)

# 绘图
plt.figure(figsize=(4, 3))
plt.plot(x_values, optimized_proof_size_F_d1*F_size + optimized_proof_size_G_d1*G_size,  marker='o')
plt.plot(x_values, optimized_proof_size_F_d2*F_size + optimized_proof_size_G_d2*G_size,  marker='o')
plt.plot(x_values, optimized_proof_size_F_d4*F_size + optimized_proof_size_G_d4*G_size,  marker='o')
plt.xlabel('m')
plt.ylabel('size (kB)')

plt.legend()
plt.grid(True)
plt.legend(frameon=False)
plt.tight_layout()


plt.show()