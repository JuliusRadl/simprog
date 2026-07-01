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

## Von oben bis unten
- Konstanten
  - alles ist in Konstanten definiert
  - auch die Beschreibung der Regimes
- **mutable struct player**
  - hat
    - Aktivierung a
    - Ein "genom" B mit:
      - Hill Saturation constant K aus [-2, -0.1] und [0.1, 2]
      - Hill cooperativity constant n zwischen [1, 9]
    - constructor
      - mini activation am Anfang
      - zufälliges K und n
- **mutable struct wattworld**
  - Definiert Koordiationsvariablen
    - Ressourcenlevel R
    - Feedrate feed
    - Zielressourcenlevel omega
    - Das feed/omega regime
    - und die zeit
  - hat private backup fields für runge kutta half steps
    - activation der spieler
    - aktuelles ressourcen level
    - aktuelle zeit
  - hat einen vector players
    - entspricht der "structure" des Agenten
  - constructor
    - performt garbage collection
    - setzt einen seed für reproducible rng
    - initialisert neue Spieler
    - Legt koordiationsvariablen fest
    - legt placeholder für die backp fields fest
- **struct snapshot**
  - zeichnet die player auf
  - und Ressourcenlevel, Feedrate, Omega und Zeit
  - benutzt dazu einen copy constructor
- ! Methoden
- **action()**
  - berechnet den absoluten Einfluss eines einzelnen Spielers
  - für die DGL zur Änderung von a
- **action(wattworld)**
  - berechnet die summe des absoluten Einflusses aller Spieler
  - auch für diese DGL
- **rk_backup!(wattworld)**
  - füllt die backup felder von dem struct wattworld mit momentanen a, R und t
- **rk_restore!**
  - stellt backup wieder her aus diesen feldern
  - benutzt pairs, um sich eine schleifenvariable aus watt.players zu schaffen
- **omega!, feed!**
  - setzen variablen omega und feed neu
  - diese setter werden nicht benutzt im ganzen modul => löschen
- **time_generator(watt)**
  - setzt omega und feed mit apply_regime!()
  - zieht zur Lesbarkeit alle a's und K's in Variablen acts und Ks
  - berechnet in activation den R-abhängigen Aktivierungsterm der DGL zur Änderung von a
  - und in inhibition den inhibitionsterm
  - benutzt dabei die hill funktion saturation(R, K)
  - berechnet die volle differenzialgleichung in da und DR
  - Könnte einzelne terme noch weiter in funktionen auslagern?
  - gibt da und DR als Tupel zurück
- **step!()**
  - performt einen RK2-Step und speichert Ergebnis gleich wieder in wattworld
  - macht backup der wattworld parameter
  - berechnet für den Halfstep Änderungsraten mit time_generator()
  - berechnet halbe schrittlánge, fügt RoC und Schrittlänge den parametern hinzu
  - berechnet dann basierend darauf die neuen Änderungsraten
  - stellt parameter aus backup wieder her
  - führt fullstep mit ganzen zeitschritt durch
- **trajectory()**
  - method header arguments stimmen nciht überein
  - nimmt eine Dauer T und entwickelt die Wattworld über diese dauer weiter
  - berechnet anzahl der Zeitschritte
  - berechnet, bei welcher zeit wir am ende landen
  - schafft step range mit konkreten zeiten, von anfang bis ende
  - intialisiert für jeden zeitschritt einen leeren snapshot
  - speichert gleich den ersten snapshot
  - iteriert dann vom zweiten bis zum letzten zeitschritt
    - variiert die wattworld bei instabilität mit vary!()
    - macht einen step!()
    - begrenzt das Ressourcenlevel auf [0, 9]
    - ersetzt alle Spieler mit a außerhalb von [0, 3] (kann man nicht gerade Fortpflanzung nennen)
    - zeichnet einen Snapshot auf
  - gibt die snapshots zurück
- **vary!()**
  - berechnet instability der wattworld (also wie weit r von omega abweicht)
  - generiert basierend darauf zufällig die Änderungen für die Genome der Player
  - benutzt die wrap funktion, um die K's im Rahmen zu halten
- **stability(watt)**
  - berechnet stabilität
  - capture radius: das limit, ab dem minimale stabilität erreicht ist (hier bei einer abweichung von 1)
  - capture concavity: exponent => wie schnell die stabilität zwischen idealwert und limit abnimmt
- **apply_regime!()**
  - setzt bei jedem step feed und omega neu, benutzt nicht die setter
- report
  - printet einen info string über die wattworld
  - und über player mit ausreichend a
- saturation()
  - Hill type saturation funktion
  - fängt Fall K == 0 ab (sollte nicht vorkommen wegen wrap())
- **demo()**
  - hat noch keine Agents eingebaut
  - plottet auch für einzelne spieler (macht keine sinn für viele spieler)
  - plottet a's für alle spieler
  - plottet K's und n's nur für spieler, die in der zweiten hälfte der simulation signifikante a's hatten (Kommentar gibt irrtümlich 1/4 an)
- **publish()**