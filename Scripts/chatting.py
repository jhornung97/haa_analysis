import numpy as np
import matplotlib.pyplot as plt

# Annahme: true_pts, true_etas und true_phis sind (4 x nEvents) Arrays für Higgs-Daughters
# Annahme: etas und phis sind (4 x nEvents) Arrays für Hadronen

def calculate_delta_r(daughter, hadron):
    delta_eta = daughter['eta'] - hadron['eta']
    delta_phi = np.abs(daughter['phi'] - hadron['phi'])
    delta_phi = np.minimum(delta_phi, 2 * np.pi - delta_phi)
    return np.sqrt(delta_eta**2 + delta_phi**2)

def has_close_hadron(daughter, reconstructed_hadrons, threshold):
    for hadron in reconstructed_hadrons:
        delta_r = calculate_delta_r(daughter, hadron)
        if delta_r < threshold:
            return True
    return False

# Annahme: true_pts, true_etas und true_phis sind (4 x nEvents) Arrays für Higgs-Daughters
# Du kannst die gewünschten Spalten entsprechend auswählen
true_pt_values = true_pts[0, :]
true_eta_values = true_etas[0, :]
true_phi_values = true_phis[0, :]

# Annahme: etas und phis sind (4 x nEvents) Arrays für Hadronen
# Zum Beispiel: reconstructed_hadrons = np.array([{'pt': ..., 'eta': ..., 'phi': ...}, ...])

threshold = 0.2
num_bins = 10  # Anzahl der Bins

# Histogramm der wahren pT-Werte erstellen
hist, bin_edges = np.histogram(true_pt_values, bins=num_bins)

# Prozentuale Anteile für jeden Bin berechnen
percentages = []
for i in range(num_bins):
    bin_indices = np.where((true_pt_values >= bin_edges[i]) & (true_pt_values < bin_edges[i + 1]))
    bin_higgs_daughters = [{'pt': true_pt_values[j], 'eta': true_eta_values[j], 'phi': true_phi_values[j]} for j in bin_indices[0]]
    
    count_with_close_hadron = 0
    total_daughters = len(bin_higgs_daughters)
    
    for daughter in bin_higgs_daughters:
        if any(calculate_delta_r(daughter, hadron) < threshold for hadron in reconstructed_hadrons):
            count_with_close_hadron += 1
    
    percentage = (count_with_close_hadron / total_daughters) * 100
    percentages.append(percentage)

# Plot erstellen
plt.plot(bin_edges[:-1], percentages, marker='o')
plt.xlabel('Wahrer Transversalimpuls (pT)')
plt.ylabel('Prozentualer Anteil mit mindestens einem Hadron innerhalb von 0.2 Delta R')
plt.title('Delta R vs. Wahrer Transversalimpuls (in Bins)')
plt.grid(True)
plt.show()
