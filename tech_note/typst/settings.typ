// Document-wide style settings — import in main.typ and call: #show: apply-styles

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

#let apply-styles(doc) = {
  set text(lang: "en")
  set text(font: "New Computer Modern")

  // Custom heading styles (no all-caps, larger sizes, better spacing)
  show heading.where(level: 1): it => block(
    above: 2em,
    below: 1em,
    text(size: 18pt, weight: "bold")[
      #if it.numbering != none {
        counter(heading).display(it.numbering)
        "."
        h(0.5em)
      }
      #it.body
    ]
  )

  show heading.where(level: 2): it => block(
    above: 1.5em,
    below: 1em,
    text(size: 16pt, weight: "bold")[
      #if it.numbering != none {
        counter(heading).display(it.numbering)
        "."
        h(0.5em)
      }
      #it.body
    ]
  )

  show heading.where(level: 3): it => block(
    above: 1.25em,
    below: 0.75em,
    text(size: 14pt, weight: "bold")[
      #if it.numbering != none {
        counter(heading).display(it.numbering)
        "."
        h(0.5em)
      }
      #it.body
    ]
  )

  doc
}
