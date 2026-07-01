# Fortpflanzung
- Mutation passiert nur hier
- Cross-Over der Genome
- Abhängig vom Energielevel und der Präsenz von Nachbarn

# Player
- haben ein Genom
- bestimmen ihr Verhalten (Ks) abhängig vom Genom und einer Variation
- Variation wird niedriger, je näher Temperatur an Soll-Wert
- Haben ein Energielevel, das kontinuierlich sinkt
- Bewegen sich weniger, wenn Temperatur besser
- bewegen sich entlang des gradienten?

# Funktionen
- 

# Agents.jl
- Welt muss zu ABM World werden mit extent
- Spieler müssen @agents werden
- brauchen kontinuierliche Plots über das observable

# Welt / Feld
- temperaturänderung muss für jede zelle berechnet werden (matrix-multiplikation)
- temperatur muss in jedem schritt auch diffundieren
- Instabilität muss lokal berechnet werden, zb in einem radius oder auch einfach pro zelle

# Interaktivität
- Statt trajectory und snapshots brauchen wir kontinuierliche entwicklung der welt
- apply_regimes können wir trotzdem behalten

# Bewertung
- brauchen sinnvolles Maß, um Erfolg der Agenten zu beurteilen, zb Anteil der Zellen mit
  der richtigen Temperatur ("habitable zones")