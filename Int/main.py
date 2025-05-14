import subprocess
import csv
from datetime import datetime
import numpy as np

# Настройки запуска
EXECUTABLE = "./main"  # путь к твоему исполняемому файлу
A = -4.99
B = 2.0
ACC = 1e-10
THREADS_LIST = np.array(2**np.linspace(0, 14, 100), dtype = np.int64)  # потоки, которые хочешь протестировать
REPEATS = 5  # сколько раз запускать для каждого числа потоков
OUTPUT_CSV = "results.csv"

def run_once(a, b, acc, threads):
    args = [EXECUTABLE, str(a), str(b), str(acc), str(threads)]
    try:
        result = subprocess.run(args, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split("\n")
        integral = float(lines[0])
        time = float(lines[1])
        return integral, time
    except subprocess.CalledProcessError as e:
        print("Ошибка при выполнении:", e.stderr)
        return None, None

def main():
    with open(OUTPUT_CSV, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["threads", "run", "integral", "time_seconds", "timestamp"])

        for threads in THREADS_LIST:
            for run in range(1, REPEATS + 1):
                integral, time = run_once(A, B, ACC, threads)
                if integral is not None:
                    writer.writerow([threads, run, integral, time, datetime.now().isoformat()])
                    print(f"[{threads} threads] Run {run}: time={time:.6f}, result={integral:.10f}")
                else:
                    print(f"[{threads} threads] Run {run} failed.")

if __name__ == "__main__":
    main()