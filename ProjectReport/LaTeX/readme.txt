
Template zur Erstellung des Projektberichts mit LaTeX an der TH Nürnberg

Kompilieren Sie das Dokument mit pdflatex => PDF

Für die Erzeugung des Literaturverzeichnisses bitte kompilieren mit 
pdflatex -> biber -> pdflatex -> pdflatex 

(alternativ processor auf bibtex umstellen, dann:
 pdflatex -> bibtex -> pdflatex -> pdflatex )


Das Template umfasst folgende Dokumente

Main.tex:
  Hauptdokument, hier soll der eigentliche Inhalt rein

Praeambel.tex
  Definitionen für das Dokumentenlayout und verwendete packages

Anhang.tex
  Dokument für sämtliche Anhänge

Title.sty
  Style-Dokument für das Titelblatt

Literatur.bib
  Database für das Literaturverzeichnis


Änderungen Vorbehalten
 Martin Michel, 18.03.2020
 Rev. 12.05.2020

