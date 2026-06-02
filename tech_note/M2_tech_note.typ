#import "@preview/texst:0.1.2": paper
#import "@preview/akatable:0.1.0": academic-table
#import "typst/settings.typ": apply-styles, btable

#show: doc => paper(
  title: [M2 Analysis],
  subtitle: [CUORE Internal Technical Note],
  authors: ((name: [Sachi Wagaarachchi]),),
  date: [#datetime.today().display("[month repr:long] [day], [year]")],
  abstract: [
    CUORE Multiplicity 2 analysis internal technial note.
  ],
  style: (
    body_font: "New Computer Modern",
    font_size: 11pt,
    page_margin: (x: 2.5cm, y: 2.5cm),
    heading_numbering: "1.1",
  ),
  doc,
)

#outline(title: [Table of Contents], indent: auto)

#pagebreak()
#outline(title: [List of Figures], target: figure.where(kind: image))

#pagebreak()
#outline(title: [List of Tables], target: figure.where(kind: table))

#show: apply-styles

#pagebreak()

#include "typst/901-prelim-notes.typ"
#pagebreak()
#include "typst/910-introduction.typ"
#pagebreak()
#include "typst/920-data.typ"
#pagebreak()
#include "typst/930-fit_model.typ"
#pagebreak()
#include "typst/940-methods.typ"
#pagebreak()
#include "typst/950-results.typ"
#pagebreak()  
#include "typst/960-systematics.typ"
#pagebreak()  
#include "typst/970-conclusion.typ"
#pagebreak()
#include "typst/990-appendix.typ"