import csv
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

CSV_FILE = "results.csv"
PLOT_FILE = "time_vs_threads.png"

def read_csv(file_path):
    times_by_threads = defaultdict(list)

    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            threads = int(row["threads"])
            time = float(row["time_seconds"])
            times_by_threads[threads].append(time)

    return times_by_threads

def plot_time_vs_threads(times_by_threads):
    threads = sorted(times_by_threads.keys())
    avg_times = [sum(times_by_threads[t]) / len(times_by_threads[t]) for t in threads]

    plt.figure(figsize=(8, 5))
    plt.plot(threads, avg_times, marker='o', linestyle='-', color='blue', label='Среднее время')
    plt.xlabel("Количество потоков")
    plt.ylabel("Время выполнения (секунды)")
    plt.title("Зависимость времени от количества потоков")
    plt.grid(True)
    plt.xticks(threads)
    plt.legend()
    plt.tight_layout()
    plt.savefig(PLOT_FILE)
    plt.show()

    a = avg_times[0]/np.array(avg_times)

    plt.figure(figsize=(8, 5))
    plt.plot(threads, a, marker='o', linestyle='-', color='blue', label='Ускорение')
    plt.xlabel("Количество потоков")
    plt.ylabel("Ускорение")
    plt.title("Зависимость Ускорения от количества потоков")
    plt.grid(True)
    plt.xticks(threads)
    plt.legend()
    plt.tight_layout()
    plt.savefig("a.png")
    plt.show()

def main():
    data = read_csv(CSV_FILE)
    plot_time_vs_threads(data)

if __name__ == "__main__":
    main()