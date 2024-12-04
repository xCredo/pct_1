import matplotlib.pyplot as plt

n_values = [10**7, 10**8]
num_threads_values = [2, 4, 6, 8]
speedup_values_n1 = [2.33, 4.6, 7, 7.8]
speedup_values_n2 = [2.45, 4.72, 6.83, 8.41]

plt.plot(num_threads_values, speedup_values_n1, label="n = 10^7")
plt.plot(num_threads_values, speedup_values_n2, label="n = 10^8")
plt.xlabel("Number of Threads")
plt.ylabel("Speedup")
plt.title("Strict Scalability Analysis")
plt.legend()
plt.savefig("mk.png")