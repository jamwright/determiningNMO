import numpy as np
import matplotlib.pyplot as plt

# loading output
data = np.genfromtxt("prob3_matter_constant.csv", delimiter=",", names=True)

LE = data["L_over_E"]
# P_ee = data["P_ee"]
# P_mumu = data["P_mumu"]
# P_mue = data["P_mue"]

# plots
channels = [
	(P_ee, r"$P(\nu_e \rightarrow \nu_e)$", "pee.png"),
	(P_mue, r"$P(\nu_\mu \rightarrow \nu_e)$", "pmue.png"),
	(P_mumu, r"$P(\nu_\mu \rightarrow \nu_\mu)$", "pmumu.png"),
	(P_mutau, r"$P(\nu_\mu \rightarrow \nu_tau)$", "pmutau.png")
]

for probability, label, filename in channels:
	plt.figure(figsize=(10,6))
	plt.plot(LE, probability, linewidth=1.5)

	plt.xlabel(r"$L/E$ (km/GeV)")
	plt.ylabel("Oscillation probability")
	plt.title(label)

	plt.ylim(0,1)
	plt.grid(True)

	plt.tight_layout()
	plt.savefig(filename, dpi=300)
	plt.show()

# combined figure plot

plt.figure(figsize=(10,6))

plt.plot(LE, P_ee, label=r"$P_{ee}$")
plt.plot(LE, P_mumu, label=r"$P_{\mu\mu}$")
plt.plot(LE, P_mue, label=r"$P_{\mu e}$")

plt.xlabel(r"$L/E$ (km/GeV")
plt.ylabel("Oscillation probability")

plt.title("Vacuum oscillations - JUNO")

plt.ylim(0,1)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("prob3_all_channels.png", dpi=300)
plt.show()

