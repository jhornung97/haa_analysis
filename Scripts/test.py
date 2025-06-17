import matplotlib.pyplot as plt

# Daten erstellen
x = range(10)
y = [i**2 for i in x]

# Plot erstellen
plt.plot(x, y)

# X-Achsen-Ticks setzen
plt.xticks(list(range(0, 10, 2)))  # Setze den ersten Tick bei 0 und dann alle zwei Einheiten

# Plot anzeigen
plt.show()
