# Meta Infos
## Nützliche Packages

- RDatasets: Port der Beispieldatensätze aus R (3500 Stück)
- Revise: Aktualisiert REPL bei Code-Änderungen (Standard bei VSC)
- Infiltrator: Debugging eigener Funktionen per @-Commands
- OhMyREPL: Syntax Highlighting im REPL
- Debugger: Einblick in fremde Funktionen per @-Commands

## Mutation Modell (Lab 114)



### Code

## Namespaces
Fehlermeldung:  "It looks like two or more modules export different bindings with this name, resulting in ambiguity";
Passiert, wenn man immer wieder include() und using .Module verwendet;
Lösung: Namespace explizit angeben

## Selection Modelling (Lab 113)

### Beobachtungen

Real ist ein Abstrakter Typ, von konkrete Typen abgeleitet
sind, die rationale Zahlen beschreiben: Int64, Float64, etc.
Aber ein Vector{Float64} ist nicht von einem Vector{Real}
abgeleitet und führt zu einem Fehler, wenn eine Methode einen
Vector{Real} fordert.

Man kann also in Methoden n::Real fordern, aber Vector{Real} 
ist unklug.

### Code

dot-Produkt:
Multipliziert a[i] mit b[i]
und bildet Summe der Produkte.
```julia
a = [1, 2, 3]
b = [2, 3, 4]
LinearAlgebra.dot(a, b)
```

## Datasets (Lab 112)

Spaltennamen eines DataFrames sehen:
```julia
names(my_df)
```

Nach Werten in einer Spalte gruppieren / sortieren:
```julia
species_groups = groupby(iris, "Species") # iris = df, "Species" = Spaltenname
```

Zeilen nach Gruppen zusammenfassen:
```julia
combine(species_groups, nrow)
```

## Chaos & Logistische Funktion (Lab 111)

Functional erstellen, also eine Funktion,
die Funktionen erstellt:
```julia
L(r) = (p -> r * p * (1-p))
```

Schauen, bei welchem Wert eine Diskrete Logistische Funktion
konvergiert:
```julia
x = 0.3
while !(L(1.1)(x) == x)
    x = L(1.1)(x)
end
```

Vektor typisieren anhand von Input:
```julia
permutations = Vector{typeof(inputVariable)}() # () anhängen, da es ja ein Konstruktor ist
```

Vektor oder Zeichen vervielfältigen:
```julia
x = [1, 2, 3]
many_x = repeat(x, 3, 3, 3) # Verfielfältige x 3-mal in der 1., 2. und 3. Dimension
```

Matrizen transposen, die nicht nur
Zahlen enthalten (erzeugt außerdem
echte Matrix, nicht nur adjoint):
```julia
x = [1 2 ; 3 4]
transposed = permutedims(x)
```

## Debugging-Tipps

Beim Debugging am Modulende Aufruf der demo()- oder unittest()-Funktion
einfügen, dann springt der Debugger beim Ausführen direkt zu Breakpoints
und man muss die Konsole nicht bemühen.

- Step Over: Zeile ausführen, ohne Funktionskörper zu betreten
- Step Into: Funktionskörper betreten
- Step Out Of: Funktionskörper verlassen
- Continue: Weitermachen bis zum Programm-Ende oder nächsten Breakpoint

Links kann man Variablen geordnet nach Lokal und Global sehen.

### Infiltrator Guide (Infiltrator.jl)

Geeignet für eigene Module, kann bei Funktionsaufrufen
nicht in fremde Pakete steigen

Funktioniert nicht immer mit dem Play Button in VSC,
eher mit include(), weil es dann im REPL läuft.

Macht Schwierigkeiten mit Makie bzw. Funktionen, die
asynchron aufgerufen werden. Deshalb zum Testen am
besten Plots rausnehmen (oder halt per include() aufrufen).

Im Code (Muss so geschrieben werden, weil man meist kein
using Infiltrator im Modul hat, und Infiltrator
nur im Main-Namensraum importiert ist):
```julia
function f(x)
	Main.Infiltrator.@infiltrate
	x^2
end

### Zeiger und Objekte

```julia
function addbody!( nbody::NBody, x0::Vector{Float64}, p0::Vector{Float64}, m::Float64=1.0)
	#= Julius' Kommentar:

	Zeile 1: Lasse nbody.x0 auf das gleiche Objekt
	wie x0 zeigen (einen Array).

	Zeile 2: Kopiere Objekt, auf das x0 zeigt und 
	lasse nbody.x darauf zeigen.

	Sinn: Wenn nbody.x später geändert wird, dann
	ändert sich nbody.x0 nicht, weil nicht beide auf
	das gleiche Objekt zeigen, auf das auch x0 beim
	Funktionsaufruf gezeigt hat. =#

	push!( nbody.x0, x0)
	push!( nbody.x, deepcopy(x0))					# Question: Why do I use deepcopy here?
	push!( nbody.p0, p0)
	push!( nbody.p, deepcopy(p0))
	push!( nbody.m, m)
	nbody.N  += 1
end
```

### Software-Entwicklung / Mocks

Entwickle Projekte von Außen nach Innen.
D. h. fange mit Ein- und Ausgabe an, also
Mock-Funktionen, die statische Rückgabe-Werte
liefern. Ein Grundgerüst aus Mock-Funktionen
ist SEHR viel Wert und ein guter Startpunkt.

Mock-Funktionen = Stubs = Fakes = Dummy Functions

## Grafik (Lab 110)

### Notizen

Man kann einen Plot nicht nachträglich in ein Axis-Objekt schieben.
Man kann aber sehr wohl ein Axis-Objekt nachträglich in ein Figure-
Objekt schieben.

Das Zeichen ";" trennt 

### Grafik

Felder eines Datentyps anzeigen lassen
(Auch wenn es mit einem Objekt des 
Datentyps nicht funktioniert):
```julia
fieldnames(Plot)
```



Attribute mit Beschreibungen einsehen
(Erst in den Help-Modus, dann nochmal
Fragezeichen):
```julia
? thing_with_attributes
```

Einfachen Graphen erzeugen:
```julia
fig = scatterlines(0:10, (0:10).^2, color=:red, linewidth=9)
```

Vergleichsfunktion für eigene Datentypen schreiben:
```julia
import Base.==
function ==( tf1::HillTF, tf2::HillTF) # Multiple Dispatch
	tf1.K==tf2.K && tf1.n==tf2.n && tf1.range==tf2.range
end
```

Axis zu Figure an bestimmter Stelle hinzufügen
```julia
fig = Figure()
ax1 = Axis(fig[1, 2])
```

Plot zu einer Axis hinzufügen:
(Entsprechend lines!(), etc)
```julia
scatterlines!(ax1, 1:10, color=:red, label="repression", linewidth=3)
```

Derzeitige Figure anzeigen:
```julia
current_figure()
```

### Observables

Observable erzeugen, Listener hinzufügen:
```julia
using Observables
t = Observable(1) # Das Observable mit dem Initialwert 1
x = map(cos, t) # Der Listener wird zum Listener, indem er ein Observable verwendet
t[] = 2
x[] # Enthält jetzt cos(2)
```

Animation durch Observable:
```julia
function animate_hill()
	tf_range = 0:30
	K = 10

	anim_n = Observable(1.0)           		# The animated cooperation value is Observable.
	hill_curve = map(anim_n) do n			# Create an expression curve that listens
		expression(HillTF(tf_range,K,n))    # directly to the current value (n) of anim_n.
	end

	fig = lines(tf_range, hill_curve)       # Create a plot of the (current) hill_curve ...
	display(fig)                        	# ... and display the result

	for n in -1:-0.1:-10                   	# Finally, animate the figure: Step through the
		anim_n[] = n                   		# values of n, and watch the hill_curve change
		sleep(0.1)                      	# as the value of anim_n moves from 1 to 10.
	end
end
```

### Clean Code (Lab 9)

Schnelles Abbruch-Kriterium mit 
kurzgeschlossenem UND-Operator:
```julia
N < 2 && return Int[]
```

### Dateisystem (Lab 8)
<!-- erzeugt einen String mit \\ (je nach OS) -->
`dataPath = joinpath(pwd(), "Subfolder1\\Subfolder2", file1.md)`

<!-- erzeugt oder öffnet Datei mit R/W-Rechten -->
`ioStream = open("ani.dat", "w")`

<!-- schreibt Text rein -->
`write(ioStream, "Mein Text")`

<!-- schließt die Datei -->
`close(ioStream)`

<!-- Hat der Zeiger das Ende der Datei erreicht? -->
`eof(ioStream)`

<!-- Daten binär auslesen -->
`data = read(ioStream)`

<!-- Binärdaten zu String konvertieren (data wird dabei geleert) -->
`str = String(data)`

<!-- Kann Hexzahlen zu Zeichen konvertieren -->
`Char(0x35)`

<!-- Datei löschen -->
`rm("ani.dat")`

<!-- Regex: ersetze 1+ Leerzeichen durch 1 Leerzeichen -->
`replace(sample_text, r"\s+" => " ")`

<!-- Splat Operator: Vektoren entpacken -->
`combinedVector = [[vektorA]..., [vektorB]...]`

<!-- String formattieren wie in Java -->
`import Printf`
`Printf.@sprintf("%30s", "20 Zeichen mit Padding)`

<!-- Formattierten String ausgeben -->
`import Printf`
`Printf.@printf("%30s", "20 Zeichen mit Padding)`

<!-- Bool Array mit True initialisieren -->
`trues(len)`

<!-- Änderungen schnell laden -->
`inc() = include("My\\Path\\MyModule.jl")`
<!-- Dann immer nach Speichern von Änderungen: -->
`inc()`
`const mm = MyModule`

### Sets / Hashing / Pairs / Tuples (Lab 7)

--- Set erzeugen
Set([1, 1, 1, 2, 3])

--- Pair abfragen
pair = 5 => 2
last(pair) => gibt 2

--- Dict erzeugen
dicty = Dict( "pi" => 3.142, "e" => 2.718)

--- Dict verändern
dicty["a"] = 23 => erzeugt ein neues Key-Value-Pair
delete!(dicty, "a")

# Doc Comments / Hilfskommentare / Documentation
```jl
"""
	Der Typ meiner Funktion

Beschreibung meiner Funktion
"""
struct myStruct # wenn ich das hover, wird der Doc Comment angezeigt
```

# Mehrzeilige Kommentare
```jl
#= Hier 
steht mein
Kommentar =#
```

# Typischer Ablauf beim Zugänglichmachen von Methoden
```jl
export myMethod
Speichern
include(<Filepath>)
using .MyModule
myMethod() # Funktioniert jetzt
```

# Zufallszahlen

```jl
r1 = rand(1:9, (2,3,4)) # 3D-Array mit ganzen Zufallszahlen zwischen 1 und 9
r2 = rand(Int, [1,2]) # Zufällig 1 oder 2
```

# Ranges
```jl
testRange = 1:1:10 # Start bei 1, Schrittlänge 1, Stop bei 10
m = reshape(1:16, 4, 4) # erzeugt 4x4 reshaped Array mit aufsteigenden Werten in Spalten != Matrix
collect(m) # macht Matrix aus reshaped Array
```

# Vektoren
```jl
# per Comprehension paarweise kombinieren
compVec2 = [x*y for x in 1:5, y in 5:10]

# Concatenation
vector1 * vector2 == prod(vector1, vector2)
size(matrix1) # Matrixdimensionen
length(matrix1) # Elementeanzahl
ndims(my_matrix) # Anzahl Dimensionen

### Vorbelegen / initialisieren
ones(3, 4)
zeros(3, 4)
fill("inhalt", 3, 4) # flexibler

--- Bestehenden Vektor füllen
b = Array{Int}(undef, 3, 4)
fill!(b, 100) => Typisierung beachten

--- Speicherplatz reservieren für Array
Array{Int}(undef, 2, 3)

--- Typ der Elemente finden
eltype(my_matrix)

--- Elemente ansteuern
my_matrix[2, 3] => Element in zweiter Zeile und 3 Spalte (da erst vertikal, dann horizontal)
my_matrix[1, end-1] => Element in erster Zeile und vorletzter Spalte
my_matrix[:, 3] => alle Elemente aus der dritten Spalte
my_matrix[2:3, end] => Elemente 2 bis 3 aus der letzten Spalte
my_matrix[begin:3, end] => Elemente 1 bis 3 aus der letzten Spalte (kinda useless)
my_matrix[5] => 5. Element, wobei erst die Zeilen runter und dann die Spalten entlang gezählt wird

--- Mit Ranges Elemente ansteuern
m2 = collect(reshape(1:81, 9, 9))
m3 = m2[1:2:end, 1:2:end] => "Die 1., 4., 7. Zeile und davon die 1., 4., 7. Spalte"

--- Form verändern
my_matrix = zeros(3, 3)
reshape(my_matrix, 1, 9) => erzeugt echte Matrix (anders als bei Ranges)

--- Elemente einer Matrix als Vektor ausgeben
m[:]

--- Spalten einer Matrix vertauschen / Permutation Index
p = [1, 3, 2, 4] => kann ich dann als Index für einen Vektor mit 4 Elementen verwenden
B = A[:, [1, 3, 2, 4]] => Indiziere alle Zeilen, sowie die Spalten 1, 3, 2, 4 => Reihenfolge ändern

--- Logical Indexing
logicalMatrix = (myMatrix .< 3) => Erzeugt "mask array", Matrix mit bool-Werten
myMatrix[logicalMatrix] => pickt nur die Werte raus, bei denen die logicalMatrix den Wert 1 hat
myMatrix.*logicalMatrix => maskiert die Matrix (erhält die Form, alle Werte, bei denen logicalMatrix 0, werden zu 0)
myMatrix[myMatrix .< 3] .= 0 => viel kompakter, schnell maskieren

### Datei einlesen

is = open("Development/Computation/pax6_hs.dat","r")

--- Zeilen lesen und Zeiger ("Cursor") bewegen:

readline(is) => String
readlines(is) => Vector{String}

--- Zeiger bewegen:
seekstart(is)

### Typisieren, Typ erzwingen, ::-Operator:

open(f::Function, args...; kwargs...)

x::Int = 5      # x muss vom Typ Int sein

### Default-Wert für Parameter angeben:

readline(io::IO=stdin; keep::Bool=false)
readline(filename::AbstractString; keep::Bool=false)

### Rumrechnen

--- Rationale Zahl
8//7
>8/7

--- Range
range(0, 2*pi, 5) => letztes Argument: Anzahl der Schritte

### Conditionals

x >= 0 ? x : -x => wie shorter if in Java

### Strings

--- Strings filtern
julia> stringVector = [1, 3, 5]
3-element Vector{Int64}:
 1
 3
 5
julia> stringTest = "ABCDEF"
"ABCDEF"
julia> stringTest[stringVector]
"ACE"

--- Alle Positionen eines Characters finden:
findall('G', dataString)

--- String-Concatenation
"Ani" * " " * "Anatta" (* = Kurzform von prod())
prod(stringVector) => Geht auch mit Vektoren, die Strings enthalten

--- Teile ersetzen
replace("xd haha", "xd" => "rofl") => Pair als zweites Argument!

--- String in Zahl verwandeln
parse(Float64, "3.14159")

--- Standardfunktionen
```jl
uppercase("xd")
lowercase("xd")
titlecase("xd") => Wird zu "Xd"
startswith("xd haha", "xd") # True
endswith("xd haha", "haha") # True
```

--- String splitten
split("xd and haha", "and") => ["xd", "haha"]

--- Substring kommt vor?
occursin("xd", "xd haha") => aus irgendeinem Grund andersrum

--- Strings mit Zwischenzeichen zusammenfuehren
join(["apples", "bananas", "pineapples"], ", ", " and ") => das letzte Argument ist das letzte Zwischenzeichen
>"apples, bananas and pineapples"

--- Variablen in Strings einfuegen (formattierte Strings)
`x = pi/2`
`println("The sine of $x is $(sin(x)).")`
>"The sine of 1.5707963267948966 is 1.0."
### Funktion auf Collection anwenden:

map(['A', 'C', 'G', 'T']) do nucleotide => nucleotide bezieht sich auf die Elemente in dem Vektor
	length(findall(nucleotide, data)) / length(data)
end

### Broadcasting / dot-Operator

log.([1, exp(1), exp(2)]) => wendet log() auf gesamte Collection
matrix1[3, :] .= 1 => Alle Elemente der dritten Zeile zu 1 machen

--- Mehrere Argumente
nin(x,y) = sin(x+y)
a = [1, 2, 3]
b = [4, 5, 6]
c = [4 5 6]
nin.(a, b) => produziert column vector, weil zwei column vector => paarweises (elementwise) Anwenden der Funktion => Hadamad-Produkt
a.*d => funkioniert auch, denn auch ein Vektor kann gebroadcastet werden
a./d
nin.(a, c) => produziert Matrix, weil row vector & column vector => outer product

### Macros

--- Zeit messen (bzw. Memory Allocations)
@time Computation.count_seq("GCAT",data_hs[1])
for i in 10:10:50
	@time fib(i)
end

--- Methoden-Quellcode einsehen
@edit sin(2*pi)

--- Welche Methode wird verwendet?
@which sin(2*pi)

### Ausdrücke / Expressions / Parsing

prog = sin(3+5/2)
parsedProg = Meta.parse(prog) => Das macht julia jedes Mal, wenn wir Code ausführen
eval(parsedProg) => Ausdruck auswerten
eval(parsedProg.args[2]) => Kann auch Teile des Ausdrucks auswerten
dump(parsedProg) => Genaue Beschreibung des Ausdrucks
parsedProg.args[2].expr => Mit . kann man Eigenschaften der Expression aufrufen

# Mehr Informationen

```jl
dump(sin) # Kann auf alles angewendet werden, liefert Nothing als Rückgabe
methods(sin) # alle Methoden einer Funktion anzeigen
fieldnames(Weasel) # "Instanzvariablen" anzeigen
supertype(Weasel) # Superklassen
subtypes(Organism) # Kindklassen
display(Weasel) # Infos über die Eigenschaften einer Variable
getfield(structname, fieldname) # hole Wert in einem Struct
```

# Konstruktoren

Abstrakten Typ definieren
<: ist ein Operator / Funktion, der auch true/false zurückgeben kann
```jl
abstract type Animal <: Organism end
```

Komplexere Typen einzeilig definieren
";" signalisiert Zeilenumbruch
```jl
struct Tree <: Organism; name::String; age::Int; end
```

# Modulmanagement

Precompilation verhindern,
Ganz oben im Modul einfüge
```jl
__precompilation__(false)
module Xyz
```

# Infos über Umgebung

Übersicht über momentane Variablen
```jl
varinfo()
```
Imports anzeigen
Zeigt mir z. B. Revise, wenn ich das in
~/.julia/config/startup.jl hinzugefügt hab
```jl
names(Main, imported=true)
```