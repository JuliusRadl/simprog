## Terminology
- Coordination variables: Umgebungsvariablen


## Im Paper
- R
  - Gemeinsame Ressource, die alle Spieler beeinflussen
  - zb Temperatur
  - ist eine "coordination variable"
  - R hängt von Aktivierung a ab
- K
  - B.K im Wattworld Module
  - Individual Saturation Attribute
  - ob ein Agent die Ressource positiv oder negativ beeinflusst
  - positive Werte führen zu Konsum und damit Reduzierung der Ressource durch den Agenten
  - Veränderbar
  - Zwischen -2 und 2, aber 0 ausgenommen
  - wird zufällig initalisiert für jeden agenten
- a
  - Aktivierung
  - Wie stark ein Agent die Ressource beeinflusst
  - zwischen 0 und 3 (0 ausgenommen)
  - a hängt von Ressource R ab
- F(t)
  - Passiver Zufluss der Ressource ins System
- beta-R
  - Natürlicher Verlust der Ressource aus dem System
  - als Koeffizient zwischen 0 und 1 eingesetzt
  - zb durch Wärmeabgabe nach außen
- x
  - Eigentlich erstes Argument der Hill-Funktion
  - Hier: Das Produkt der Aktivierung a und dem Betrag des Saturation level K
    summiert über alle Agenten
  - Also je mehr Agenten aktiv sind, desto mehr wird jeder einzelne Agent dadurch
  beeinflusst
- K'
- h
  - Hill Saturation Funktion
  - Bei negativem K nehmen wir einfach den Betrag von K und ziehen dann h(x, K) 
  von 1 ab => Aktivierung wird sehr klein
  - 
- delta Kp
  - Zufällige Variation der individuellen Saturation K
  - nach jedem Zeitschritt
  - Maximal so hoch, wie die Abweichung zwischen R und Omega, aber nur bis zu einer
  Abweichung von 1 (nimmt Minimum)
  - mal 0.05
  - Dh bei großer Abweichung ist eine geringe Änderung möglich, bei geringer Abweichung
  eine sehr geringe Änderung und bei keiner Abweichung gar keine Änderung
  - Minimal 0
- Omega
  - Zielwert der Ressource
  - zb Idealtemperatur
- n
  - Im Paper "Cooperativity" genannt
  - in Wattworld B.n
  - wird in WW zufällig initialisiert
  - wird auch nach jedem Zeitschritt zufällig variiert
  - Exponent der drei Terme in der Hillfunktion
  - Bestimmt die Steilheit der Hillfunktion
  - Bei hoher Cooperativity bleibt die Aktivierung des Agenten lange niedrig und
  steigt dann rapide an, wenn
- beta-a
  - Natürlicher Verlust der Aktivierung jedes Agenten

## Feedback Palfreyman
- Daisyworld Paper lesen
- Temperatur ist lokal
- Agent hat keine Absicht, versucht nicht, irgendwas zu erreichen.
- Wir geben nur das Stabilisierungskriterium vor: Unter der und der örtlichen Feldausprägung
  sinkt die persönliche Variation des Agenten auf Null (also wenn keine Abweichung mehr da ist)
- Wenn Struktur das Feld bestimmt und umgekehrt, dann ist das Stabilisierung. Die Fragen ist,
  ob diese Stabilisierung zu Adaption oder Modularität etc. führt. Dann haben wir lebensähnliche
  Zustände gezeigt.
- Die Frage ist, was dann das Leben ist: 
- Tipp: Auf jeden Fall eine Störung irgendeiner Art einbauen, um Informationen über die
  Fähigkeit zur Adaption zu erhalten (interessant und einfach): Bleiben bestehende "Leben" erhalten
  und adaptieren ihr Verhalten oder entstehen neue? (Einfärben!)
- Braucht ausreichend Agenten nah beieinander, damit sie sich gegenseitig beeinflussen können
- Bei Bewegung: Darf nicht vorgeben, dass Agenten sich zb Richtung steigenden Gradienten bewegen:
  Agenten müssen dieses Verhalten selbst erlernen

## Plan
- Lokales Temperaturfeld
- Temperatur diffundiert und nimmt natürlich ab
- Agenten bewegen sich
- Agenten konsumieren und erzeugen die Ressource