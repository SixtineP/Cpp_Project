import csv
import matplotlib.pyplot as plt

CSV_FILE = "Convergence_results.csv"

def read_csv(path):
    data = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({
                "n": int(row["n"]),
                "m": int(row["m"]),
                "err": float(row["err"])
            })
    return data

def plot_convergence(data, ratio=2):
    pts = [(d["n"], d["err"]) for d in data if d["m"] == ratio * d["n"]]
    pts.sort(key=lambda x: x[0])

    n_vals = [p[0] for p in pts]
    err_vals = [p[1] for p in pts]

    plt.figure(figsize=(9, 6))
    plt.plot(n_vals, err_vals, marker="+")
    plt.xlabel("Number of time steps n (m = 2n)")
    plt.ylabel("Error |PDE - BS|")
    plt.title("Convergence of PDE price to Black–Scholes price")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("Convergence.png")
    plt.close()

def main():
    data = read_csv(CSV_FILE)
    plot_convergence(data, ratio=2)

if __name__ == "__main__":
    main()
