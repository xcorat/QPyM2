#import "@preview/texst:0.1.2": paper
#import "@preview/akatable:0.1.0": academic-table
#set text(lang: "en")

#show: doc => paper(
  title: [M2 Analysis],
  authors: ((name: [Sachi Wagaarachchi]),),
  date: [#datetime.today().display("[month repr:long] [day], [year]")],
  abstract: [
    CUORE Multiplicity 2 analysis note.
  ],
  style: (
    body_font: "New Computer Modern",
    font_size: 11pt,
    page_margin: (x: 2.5cm, y: 2.5cm),
    heading_numbering: "1.1",
  ),
  doc,
)

#outline(title: [List of Figures], target: figure.where(kind: image))
#outline(title: [List of Tables], target: figure.where(kind: table))

#set text(font: "New Computer Modern")

#let btable(..args) = block(stroke: (bottom: 1.5pt))[
  #set table(
    stroke: (x, y) => (
      left: none,
      right: none,
      top: if y == 0 { 1.5pt } else if y == 1 { 0.75pt } else { none },
      bottom: none,
    ),
    inset: (x: 8pt, y: 5pt),
    align: left,
  )
  #table(..args)
]

#let ndbd = "ndbd" //[$0 nu beta beta$]

#pagebreak()

#include "typst/910-introduction.typ"
#pagebreak()
#include "typst/920-data.typ"
#pagebreak()
#include "typst/930-fit_model.typ"
#pagebreak()
#include "typst/940-model_selection.typ"
#pagebreak()
#include "typst/950-results.typ"
#pagebreak()  
#include "typst/960-systematics.typ"
#pagebreak()
#include "typst/990-appendix.typ"