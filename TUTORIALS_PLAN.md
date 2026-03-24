# Plan zajęć: Spatiotemporalna propagacja sygnału

**Biologia Systemów · Uniwersytet Warszawski**

Paczka zawiera skrypty analityczne, notatki matematyczne (handout + arkusz ćwiczeniowy w formacie PDF) oraz metadane eksperymentu. Zajęcia podzielone są na trzy tutoriale: pierwszy poświęcony jest eksploracji danych, drugi — samodzielnej pracy ze skryptami i interpretacji wyników statystycznych, trzeci — krytycznej analizie ograniczeń metody i implementacji własnych ulepszeń. Dane do analizy (ścieżki pojedynczych komórek, pliki `.csv.gz`) należy pobrać osobno z sekcji *Dane* na Moodlu i umieścić w głównym katalogu projektu. Przed pierwszymi zajęciami warto zapoznać się z handoutem matematycznym (`notes/spatiotemporal_propagation_handout.pdf`), który wprowadza całą notację używaną w skryptach i arkuszach.

---

## Tutorial 1 — Eksploracja danych

*Folder: `notebooks-part1/`*

Pięć notebooków wprowadzających do datasetu. Studenci pobierają i wstępnie przetwarzają dane, wizualizują dynamikę sygnałów ERK i AKT w pojedynczych komórkach, analizują wzorce czasowe w ścieżkach obserwacji, badają przestrzenną organizację komórek w polach widzenia oraz identyfikują przejawy kolektywnej aktywności w monowarstwach.

| Notebook | Temat |
|---|---|
| `00_get_data.ipynb` | Dostęp do danych i wstępne przetwarzanie |
| `01_explore_erk_akt_signaling.ipynb` | Rozkłady sygnałów ERK i AKT |
| `02_temporal_dynamics.ipynb` | Dynamika sygnału w czasie wzdłuż ścieżek |
| `03_spatial_structure.ipynb` | Przestrzenna organizacja komórek |
| `04_collective_behavior.ipynb` | Kolektywne wzorce aktywności |

---

## Tutorial 2 — Analiza propagacji

*Folder: `notebooks-part2/` · Part I arkusza ćwiczeniowego*

Trzy notebooki odpowiadające blokom Part I arkusza. Studenci uruchamiają skrypt analityczny, interpretują wyjściowe statystyki i badają wpływ parametrów na wyniki.

| Notebook | Zadania | Zakres |
|---|---|---|
| `Part1_Block1_FirstRun.ipynb` | 1–2 | Uruchomienie skryptu; odczyt i weryfikacja statystyk (węzły, krawędzie, próg θ, RR) |
| `Part1_Block2_ParameterSensitivity.ipynb` | 3–5 | Wrażliwość na promień *r*, okno czasowe *W* i kwantyl progu *q* |
| `Part1_Block3_Comparison.ipynb` | 6–7 | Porównanie dwóch warunków (WT vs. mutant); analiza kohortowa wszystkich mutacji |

---

## Tutorial 3 — Ulepszenia metody

*Folder: `notebooks-part3/` · Part II arkusza ćwiczeniowego*

Trzy notebooki odpowiadające ograniczeniom L1–L3 opisanym w arkuszu. Studenci identyfikują matematyczne słabości bazowej metody, wyprowadzają ulepszone estymatory i implementują je w Pythonie.

| Notebook | Ograniczenie | Ulepszenie |
|---|---|---|
| `Part2_L1_LaggedExposure.ipynb` | L1 — simultaniczność ekspozycji i odpowiedzi | Opóźniona ekspozycja $E_k^{(\tau)}(t)$; krzywa $\mathrm{RR}^{(\tau)}$ w funkcji opóźnienia $\tau$ |
| `Part2_L2_GlobalThreshold.ipynb` | L2 — nieporównywalność progów między blokami | Globalny próg $\theta^*$ obliczony ze wszystkich bloków; porównanie wyników przed i po korekcie |
| `Part2_L3_DoseResponse.ipynb` | L3 — binarna ekspozycja pomija intensywność sąsiedztwa | Ułamkowa ekspozycja $\tilde{E}_k(t)$; krzywa dawka–odpowiedź i gradowane ryzyko względne $\widetilde{\mathrm{RR}}$ |

---

## Materiały do przeczytania przed zajęciami

| Plik | Opis |
|---|---|
| `notes/spatiotemporal_propagation_handout.pdf` | Formalizacja matematyczna potoku analitycznego: grafy spatiotemporalne, detekcja skoków, estymatory ekspozycji i propagacji, miary ryzyka |
| `notes/classroom_worksheet.pdf` | Arkusz ćwiczeniowy na 3 godziny: zadania prowadzone (Part I) i otwarte zadania implementacyjne (Part II) |
