# Transcription of Huynh's (1977) data from Choi (2013)

Source: `reading/Choi2013.pdf`, **printed page 25** (= PDF page 32), Figure 1,
captioned "Figure 1. Huynh's (1977) data". The surrounding text reads: "A toy
data set (Huynh, 1977) shown in Figure 1 is used for this study. This synthetic
data contains 12 persons and 6 items with random missing responses."

## Method

The page was rendered to PNG at **300 dpi** with

```
pdftoppm -f 32 -l 32 -r 300 -png reading/Choi2013.pdf choi_p25
```

and the printed page number "25" at the foot of the rendering was checked
against the printed page cited in the revision plan. The matrix block was then
re-rendered at **600 dpi** over the figure's bounding box,

```
pdftoppm -f 32 -l 32 -r 600 -png -x 1900 -y 580 -W 1400 -H 1450 \
         reading/Choi2013.pdf choi_matrix_zoom
```

and every cell was read off that enlargement. The figure is set in a monospaced
face, one row per person and six columns, with missing responses printed as
`NA`.

## Cells that could not be read

**None.** All 72 cells (12 persons x 6 items) were unambiguous at 300 dpi and
again at 600 dpi. There is no smudged, clipped or ambiguous entry, and no cell
required a judgement call between `0`, `1` and `NA`.

## Independent cross-check

The PDF carries a text layer for this figure. `pdftotext -f 32 -l 32` returns the
twelve rows verbatim, and they agree with the visual transcription cell for cell.
The two readings are independent (glyph rendering versus embedded character
codes) and neither was used to correct the other.

## The transcribed matrix

```
person   i1  i2  i3  i4  i5  i6
  1       1   1   1   1   0   1
  2       1   0   1   1   1   1
  3      NA   1   1   1   0   0
  4       0  NA  NA   1   0   1
  5       1  NA   0   1   0  NA
  6       1   1   1   1   0   1
  7       0  NA   1  NA   0   0
  8       0  NA   1   1   1   1
  9       1   1   0   1  NA   0
 10      NA   0   0   0  NA   0
 11       1   1   1   1   1   1
 12      NA   0  NA   0   0   1
```

13 of the 72 cells are missing, leaving **59 observed responses**. This is the
unbalanced, incomplete p x i design Brennan (2001, p. 227) uses.

Recorded as `huynh1977.csv` in wide form, one row per person, `NA` for missing.
