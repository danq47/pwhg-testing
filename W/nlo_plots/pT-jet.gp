# ======================================================
# - Plot 2/6: pT-jet 
# ======================================================

reset

# GnuPlot index: 
# ---------------
theIndex = 1

# Input files: 
# -------------
file_1  = "NLO.top"
file_2  = "LHEF.top"

# Setting colours, line-types and pointer types: 
# -----------------------------------------------
colour_1  = "black"
colour_2  = "red"
lt_1      = "1"
lt_2      = "2"
pt_1      = "0"
pt_2      = "2"

# Titles for legend (generated from input files): 
# ------------------------------------------------
title_1 = "NLO"
title_2 = "LHEF"

# Titles for x-axis and y-axis: 
# ------------------------------
thexlabel = "pT-jet [units]"
set label \
'${\mathrm d}\sigma/{\mathrm d}{\mathrm O}$ [pb/units]' \
at screen 0.10 , 0.45 rotate
set label "Ratio" at screen 0.10 , 0.22 rotate

# Settings for the legend: 
# -------------------------
# set key left / right / top / bottom / outside / inside
set key right top
set key spacing 1.6

# Data combined pairwise from input files: 
# -----------------------------------------
both_1  = '< paste NLO.top NLO.top'
both_2  = '< paste NLO.top LHEF.top'

# Output file name based on title: 
# ---------------------------------
set output 'pT-jet.tex'

# Blurb to set up axes, margins etc: 
# -----------------------------------
set terminal epslatex color
set multiplot

set origin 0.0,0.0
set size   1.0,1.0
set lmargin at screen 0.25
set rmargin at screen 0.95
set tmargin at screen 0.95
set bmargin at screen 0.40

set datafile fortran
# Plot as lines or points?
set style data step
# set style data xyerrorbars

set xtics format ""
set xrange [0:300.00000000] writeback

set macros

# To plot points use: 
# --------------------
range1 = "(($1+$2)/2):3:1:2:($3-$4):($3+$4)"
range2 = "($1+$2)/2:($7/$3):1:2:($7/$3)*(1-(($4/$3)**2+($8/$7)**2)**0.5):($7/$3)*(1+(($4/$3)**2+($8/$7)**2)**0.5)" 

# To plot lines use: 
# -------------------
range1L = "(($1+$2)/2):3" 
range2L = "(($1+$2)/2):($7/$3)" 

set log y
unset xtics

# The bit that does the main plots: 
# ----------------------------------
plot \
file_1 index theIndex using @range1L \
title title_1  \
lc rgb colour_1 lw 1 lt lt_1 \
,\
file_2 index theIndex using @range1L \
title title_2  \
lc rgb colour_2 lw 1 lt lt_2 \


# Blurb to set up the ratio plot window: 
# ---------------------------------------
set tmargin at screen 0.4
set bmargin at screen 0.2
set nolog y
set yrange [0.:2]
set xrange restore
set xlabel thexlabel
set format x
set key off
set ytics("0" 0,"0.5" 0.5,"1" 1,"1.5" 1.5," " 2)
set xtics

# The bit that does the ratio plots: 
# -----------------------------------
plot \
both_2 index theIndex using @range2L title "" \
lc rgb colour_2 lw 3 lt lt_2 \
,\
1 with lines lt 2 lc 0


# The last bit: 
# --------------
unset multiplot

shell0 = "sed -i -e \"s/input{.*}/input{pT-jet.tex}/\" tot.tex "
shell1 = "dvipdf tot.dvi pT-jet.pdf"
shell2 = "pdfcrop pT-jet.pdf"
shell3 = "dvips -E tot.dvi -o tot.eps"
shell4 = "mv pT-jet-crop.pdf pT-jet.pdf"
shell5 = "mv tot.eps pT-jet.eps"
system shell0
set output
system "latex tot.tex"
system shell1
system shell2
system shell3
system shell4
system shell5

