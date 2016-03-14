#############################################
#!/bin/bash                          # blah #
theOS=`uname`                        # blah #
today=`date "+%d_%m_%y"`             # blah #
if [ "$theOS" = "Darwin" ]           # blah #
then thesed=sed                     # blah #
else thesed=sed                      # blah #
fi                                   # blah #
#############################################
sed -i '' 's/_/-/g' *.top;
# 0. Set a name for the gnuplot script which we will generate and,
#    more usefully, a prefix for the eps files for the plots we will
#    be generating. Also you can set whether you want "eps" or "pdf"
#    output files with the outputType variable..
prefix=""         # <--- prefix for all produced .gp files
outputType="eps"  # <--- eps for eps and pdf or pdf for just pdf
plotPoints="n"    # <--- y for points, n for lines

echo
echo "*********************"
echo "**** genplots.sh ****"
echo "*********************"

#
# 1. Check we got input at least two data files:
#
if [ $# -lt 2 ]
then
    echo
    echo "Usage: genplot.sh file1.top ... fileX.top (X<10)."
    echo "Exiting ..." 
    echo
    exit
else
    nInputFiles=$#
    echo
    echo "Found $nInputFiles input files on the command line."
    echo "-----------------------------------------"
    echo ""
fi

# <----- 9 is the max number of files ...

#
# 2. Store the input file names:
#
file=(dummy $1 $2 $3 $4 $5 $6 $7 $8 $9)
for (( ixx=1 ; ixx<=$nInputFiles ; ixx++ ))
do
    echo "${file[ixx]}"
done
echo

#
# 3. Some default colour, line-type and pointer-type settings
#
colour=(none black red blue dark-green purple turquoise magenta cyan yellow orange violet)
#colour=(none red red red dark-green dark-green dark-green blue blue blue)
pt=(0 0 2 3 4 5 6 7 8 9)
lt=(0 1 2 4 1 2 4 1 2 4)

#
# 3. Loop over the histogram titles in file[1] -- tagged by
#    #'s -- and save these in the array "titles":  
#
ixx=0
for i in `$thesed -n "/#/ =" ${file[1]}`
do
  ixx=$((ixx+1))
  titles[$ixx]=`$thesed -n "$i p" ${file[1]} | $thesed "s/#//"`
  titles[$ixx]=`echo ${titles[ixx]} | sed "s/ index \(.*\)//"`
  theLineNo=`$thesed -n "$i =" ${file[1]}`
  theLineNo=$((theLineNo+1))
  theLine=`$thesed -n "$theLineNo p" ${file[1]}`
  theLine=`echo $theLine | $thesed "s/\(.*\) \(.*\) \(.*\) \(.*\)/\1/"`
  theLine=`echo $theLine | $thesed "s/D+\([0-9]\)\([0-9]\)/*10^(\1\2)/"`
  theLine=`echo $theLine | $thesed "s/D-\([0-9]\)\([0-9]\)/\/10^(\1\2)/"`
  theLine=`echo $theLine | $thesed "s/E+\([0-9]\)\([0-9]\)/*10^(\1\2)/"`
  theLine=`echo $theLine | $thesed "s/E-\([0-9]\)\([0-9]\)/\/10^(\1\2)/"`
  theLine=`echo "scale=5 ; $theLine " | bc`
  xMin[$ixx]=$theLine
done
ixx=0
jxx=0
for i in `$thesed -n "/^$/ =" ${file[1]}`
do
  ixx=$((ixx+1))
  if [ $((ixx%2)) = 0 ] ; then continue ; fi
  jxx=$((jxx+1))
  theLineNo=`$thesed -n "$i =" ${file[1]}`
  theLineNo=$((theLineNo-1))
  theLine=`$thesed -n "$theLineNo p" ${file[1]}`
  theLine=`echo $theLine | $thesed "s/\(.*\) \(.*\) \(.*\) \(.*\)/\2/"`
  theLine=`echo $theLine | $thesed "s/D+\([0-9]\)\([0-9]\)/*10^(\1\2)/"`
  theLine=`echo $theLine | $thesed "s/D-\([0-9]\)\([0-9]\)/\/10^(\1\2)/"`
  theLine=`echo $theLine | $thesed "s/E+\([0-9]\)\([0-9]\)/*10^(\1\2)/"`
  theLine=`echo $theLine | $thesed "s/E-\([0-9]\)\([0-9]\)/\/10^(\1\2)/"`
  theLine=`echo "scale=5 ; $theLine " | bc`
  xMax[$jxx]=$theLine
done
#
# 4. Work out the output file names and x-ranges for each
#    of the histograms in "titles":
#
for (( ixx=1 ; ixx<=${#titles[@]} ; ixx++ ))
do
  outputFile[$ixx]=`echo ${titles[ixx]}      | $thesed "s/\.//g"   `
  outputFile[$ixx]=`echo ${outputFile[$ixx]} | $thesed "s/\=//g"   `
  outputFile[$ixx]=`echo ${outputFile[$ixx]} | $thesed "s/|//g"    `
  outputFile[$ixx]=`echo ${outputFile[$ixx]} | $thesed "s/<//g"    `
  outputFile[$ixx]=`echo ${outputFile[$ixx]} | $thesed "s/\,//g"   `
  outputFile[$ixx]=`echo ${outputFile[$ixx]} | $thesed "s/\ /\_/g" `
  if [ "$prefix" = "" ] ; then
    outputFile[$ixx]=`echo ${outputFile[$ixx]} | $thesed "s/\ /\_/g" `
  else
    outputFile[$ixx]=`echo $prefix\_${outputFile[$ixx]} | $thesed "s/\ /\_/g" `
  fi
  xRange[$ixx]=[${xMin[ixx]}:${xMax[ixx]}]
done

#
# 5. Now reassure (hopefully) with informative output:
#    
echo
echo "Found ${#titles[@]} plots in file ${file[1]}: "
echo "----------------------------------------------"
echo ""
for (( ixx=1 ; ixx<=${#outputFile[@]} ; ixx++ ))
do
  echo $ixx":" ${titles[ixx]} ", range " ${xMin[ixx]} "->" ${xMax[ixx]}
done
echo 

#
# 7. Create the file needed for getting LaTeX output:
#
echo "\documentclass[14pt]{extarticle}% one can choose 8,9,10,11,12,14,17,20 pt" > tot.tex
echo "\usepackage{graphicx,color}" >> tot.tex
echo "\begin{document}" >> tot.tex
echo "\thispagestyle{empty}" >> tot.tex
echo "\input{data.tex}" >> tot.tex
echo "\end{document}" >> tot.tex

#
# 8. Loop over all the different plots:
#
echo
echo "Writing GnuPlot script: "
echo "------------------------"
echo ""
for (( ixx=1 ; ixx<=${#titles[@]} ; ixx++ ))
do
#
# 6. Create the file for the output gnuplot script:
#
  if [ "$prefix" = "" ] ; then
      gnuPlotScript=${titles[ixx]}.gp
  else
      gnuPlotScript=$prefix\_${titles[ixx]}.gp
  fi
  if [ -e $gnuPlotScript ]
  then
      rm $gnuPlotScript
  fi
  touch $gnuPlotScript
  if [ $((ixx%10)) -eq 0 ] ; then echo "Written script for $ixx/${#titles[@]} plots ..." ; fi
#
# 9. Write the blurb for the ixx-th plot:
#
  echo "# ======================================================" >> $gnuPlotScript 
  echo "# - Plot $ixx/${#titles[@]}: ${titles[ixx] } "            >> $gnuPlotScript
  echo "# ======================================================" >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "reset"                                                    >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "# GnuPlot index: "                                        >> $gnuPlotScript
  echo "# ---------------"                                        >> $gnuPlotScript
  echo "theIndex = $((ixx-1))"                                    >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "# Input files: "                                          >> $gnuPlotScript
  echo "# -------------"                                          >> $gnuPlotScript
  for (( jxx=1 ; jxx<=$nInputFiles ; jxx++ ))
  do
    echo "file_$jxx  = \"${file[jxx]}\""                          >> $gnuPlotScript
  done
  echo ""                                                         >> $gnuPlotScript
  echo "# Setting colours, line-types and pointer types: "        >> $gnuPlotScript
  echo "# -----------------------------------------------"        >> $gnuPlotScript
  for (( jxx=1 ; jxx<=$nInputFiles ; jxx++ ))
  do
    echo "colour_$jxx  = \"${colour[jxx]}\""                      >> $gnuPlotScript
  done
  for (( jxx=1 ; jxx<=$nInputFiles ; jxx++ ))
  do
    echo "lt_$jxx      = \"${lt[jxx]}\""                          >> $gnuPlotScript
  done
  for (( jxx=1 ; jxx<=$nInputFiles ; jxx++ ))
  do
    echo "pt_$jxx      = \"${pt[jxx]}\""                          >> $gnuPlotScript
  done
  echo ""                                                         >> $gnuPlotScript
  echo "# Titles for legend (generated from input files): "       >> $gnuPlotScript
  echo "# ------------------------------------------------"       >> $gnuPlotScript
  for (( jxx=1 ; jxx<=$nInputFiles ; jxx++ ))
  do
    humanReadableTitle=`echo ${file[jxx]} | $thesed "s/_/ /g" | $thesed "s/.top//g"`
    echo "title_$jxx = \"$humanReadableTitle\""                   >> $gnuPlotScript
#   echo "title_$jxx = \"${file[jxx]}\""                          >> $gnuPlotScript
  done
  echo ""                                                         >> $gnuPlotScript
  echo "# Titles for x-axis and y-axis: "                         >> $gnuPlotScript
  echo "# ------------------------------"                         >> $gnuPlotScript
  echo "thexlabel = \"${titles[ixx]} [units]\""                   >> $gnuPlotScript
  echo "set label \\"                                             >> $gnuPlotScript
  echo "'\${\\mathrm d}\\sigma/{\\mathrm d}{\\mathrm O}\$ [pb/units]' \\" >> $gnuPlotScript
  echo "at screen 0.10 , 0.45 rotate"                             >> $gnuPlotScript
  echo "set label \"Ratio\" at screen 0.10 , 0.22 rotate" >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "# Settings for the legend: "                              >> $gnuPlotScript
  echo "# -------------------------"                              >> $gnuPlotScript
  echo "# set key left / right / top / bottom / outside / inside" >> $gnuPlotScript
  echo "set key right top"                                        >> $gnuPlotScript
  echo "set key spacing 1.6"                                      >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "# Data combined pairwise from input files: "              >> $gnuPlotScript
  echo "# -----------------------------------------"              >> $gnuPlotScript
  for (( jxx=1 ; jxx<=$nInputFiles ; jxx++ ))
  do
    echo "both_$jxx  = '< paste ${file[1]} ${file[jxx]}'"         >> $gnuPlotScript
  done
  echo ""                                                         >> $gnuPlotScript 
  echo "# Output file name based on title: "                      >> $gnuPlotScript 
  echo "# ---------------------------------"                      >> $gnuPlotScript
  echo "set output '${outputFile[ixx]}.tex'"                      >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript 
  echo "# Blurb to set up axes, margins etc: "                    >> $gnuPlotScript
  echo "# -----------------------------------"                    >> $gnuPlotScript
  echo "set terminal epslatex color"                              >> $gnuPlotScript
  echo "set multiplot"                                            >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "set origin 0.0,0.0"                                       >> $gnuPlotScript
  echo "set size   1.0,1.0"                                       >> $gnuPlotScript
  echo "set lmargin at screen 0.25"                               >> $gnuPlotScript
  echo "set rmargin at screen 0.95"                               >> $gnuPlotScript
  echo "set tmargin at screen 0.95"                               >> $gnuPlotScript
  echo "set bmargin at screen 0.40"                               >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "set datafile fortran"                                     >> $gnuPlotScript
  echo "# Plot as lines or points?"                               >> $gnuPlotScript
  if [ "plotPoints" = "y" ] ; then
    echo "# set style data lines"                                 >> $gnuPlotScript
    echo "set style data xyerrorbars"                             >> $gnuPlotScript
  else
    echo "set style data step"                                   >> $gnuPlotScript
    echo "# set style data xyerrorbars"                           >> $gnuPlotScript
  fi
  echo ""                                                         >> $gnuPlotScript
  echo "set xtics format \"\""                                    >> $gnuPlotScript
  echo "set xrange ${xRange[ixx]} writeback"                      >> $gnuPlotScript
#  echo "set format y \"\$10^{%L}$\""                              >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "set macros"                                               >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "# To plot points use: "                                   >> $gnuPlotScript
  echo "# --------------------"                                   >> $gnuPlotScript
  echo "range1 = \"((\$1+\$2)/2):3:1:2:(\$3-\$4):(\$3+\$4)\""     >> $gnuPlotScript
  echo "range2 = \"(\$1+\$2)/2:(\$7/\$3):1:2:(\$7/\$3)*(1-((\$4/\$3)**2+(\$8/\$7)**2)**0.5):(\$7/\$3)*(1+((\$4/\$3)**2+(\$8/\$7)**2)**0.5)\" "                    >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "# To plot lines use: "                                    >> $gnuPlotScript
  echo "# -------------------"                                    >> $gnuPlotScript
  echo "range1L = \"((\$1+\$2)/2):3\" "                           >> $gnuPlotScript
  echo "range2L = \"((\$1+\$2)/2):(\$7/\$3)\" "                   >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "set nolog y"                                                >> $gnuPlotScript
  echo "unset xtics"                                              >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "# The bit that does the main plots: "                     >> $gnuPlotScript
  echo "# ----------------------------------"                     >> $gnuPlotScript
  echo "plot \\"                                                  >> $gnuPlotScript
  for (( jxx=1 ; jxx<=$nInputFiles ; jxx++ ))
  do
    if [ "$plotPoints" = "y" ] ; then
      echo "file_$jxx index theIndex using @range1 \\" >> $gnuPlotScript
    else
      echo "file_$jxx index theIndex using @range1L \\" >> $gnuPlotScript
    fi
    echo "title title_$jxx  \\"                                   >> $gnuPlotScript
    if [ $jxx -eq 1 ]
    then
      if [ "$plotPoints" = "y" ] ; then
        echo "lc rgb colour_$jxx lw 1 lt lt_$jxx pt pt_$jxx \\"   >> $gnuPlotScript
      else
        echo "lc rgb colour_$jxx lw 1 lt lt_$jxx \\"              >> $gnuPlotScript
      fi
    else
      if [ "$plotPoints" = "y" ] ; then
        echo "lc rgb colour_$jxx lw 1 lt lt_$jxx pt pt_$jxx \\"   >> $gnuPlotScript
      else
        echo "lc rgb colour_$jxx lw 1 lt lt_$jxx \\"              >> $gnuPlotScript
      fi
    fi
    if [ $jxx -ne $nInputFiles ]
    then
      echo ",\\"                                                  >> $gnuPlotScript
    fi
  done
  echo ""                                                         >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "# Blurb to set up the ratio plot window: "                >> $gnuPlotScript
  echo "# ---------------------------------------"                >> $gnuPlotScript
  echo "set tmargin at screen 0.4"                                >> $gnuPlotScript
  echo "set bmargin at screen 0.2"                                >> $gnuPlotScript
  echo "set nolog y"                                              >> $gnuPlotScript
  echo "set yrange [0.95:1.05]"                                         >> $gnuPlotScript
  echo "set xrange restore"                                       >> $gnuPlotScript
  echo "set xlabel thexlabel"                                     >> $gnuPlotScript
  echo "set format x"                                             >> $gnuPlotScript
  echo "set key off"                                              >> $gnuPlotScript
  echo "set ytics(\"0\" 0,\"0.5\" 0.5,\"1\" 1,\"1.5\" 1.5,\" \" 2)" >> $gnuPlotScript
  echo "set xtics"                                                >> $gnuPlotScript
#  echo "set arrow from ${xMin[ixx]},1 to ${xMax[ixx]},1 nohead lt -1 lw 0.5" >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "# The bit that does the ratio plots: "                    >> $gnuPlotScript
  echo "# -----------------------------------"                    >> $gnuPlotScript
  echo "plot \\"                                                  >> $gnuPlotScript
  for (( jxx=2 ; jxx<=$nInputFiles ; jxx++ ))
  do
    if [ "$plotPoints" = "y" ] ; then
      echo "both_$jxx index theIndex using @range2 title \"\" \\"   >> $gnuPlotScript
    else
      echo "both_$jxx index theIndex using @range2L title \"\" \\"  >> $gnuPlotScript
    fi
    if [ "$plotPoints" = "y" ] ; then
      echo "lc rgb colour_$jxx lw 3 lt lt_$jxx pt pt_$jxx \\"     >> $gnuPlotScript
    else
      echo "lc rgb colour_$jxx lw 3 lt lt_$jxx \\"                >> $gnuPlotScript
    fi
    if [ $jxx -ne $nInputFiles ]
    then
      echo ",\\"                                                  >> $gnuPlotScript
    else
      echo ",\\"                                                  >> $gnuPlotScript
      echo "1 with lines lt 2 lc 0"                               >> $gnuPlotScript
    fi
  done
  echo ""                                                         >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "# The last bit: "                                         >> $gnuPlotScript
  echo "# --------------"                                         >> $gnuPlotScript
  echo "unset multiplot"                                          >> $gnuPlotScript
  echo ""                                                         >> $gnuPlotScript
  echo "shell0 = \"$thesed -i -e \\\"s/input{.*}/input{${outputFile[ixx]}.tex}/\\\" tot.tex \""  >> $gnuPlotScript
  echo "shell1 = \"dvipdf tot.dvi ${outputFile[ixx]}.pdf\""       >> $gnuPlotScript
  echo "shell2 = \"pdfcrop ${outputFile[ixx]}.pdf\""              >> $gnuPlotScript
  echo "shell3 = \"dvips -E tot.dvi -o tot.eps\""                 >> $gnuPlotScript
  echo "shell4 = \"mv ${outputFile[ixx]}-crop.pdf ${outputFile[ixx]}.pdf\"" >> $gnuPlotScript
  echo "shell5 = \"mv tot.eps ${outputFile[ixx]}.eps\""           >> $gnuPlotScript
  echo "system shell0"                                            >> $gnuPlotScript
  echo "set output"                                               >> $gnuPlotScript
  echo "system \"latex tot.tex\""                                 >> $gnuPlotScript
  echo "system shell1"                                            >> $gnuPlotScript
  echo "system shell2"                                            >> $gnuPlotScript
  if [ "$outputType" = "eps" ]
  then
      echo "system shell3"                                        >> $gnuPlotScript
      echo "system shell4"                                        >> $gnuPlotScript
      echo "system shell5"                                        >> $gnuPlotScript
  fi
  echo ""                                                         >> $gnuPlotScript

done

echo ""
echo ""
echo "newgenplots.sh done!"
echo "--------------------"
echo ""
echo ""
echo "To generate .eps / .pdf output now type: gnuplot <file>.gp "
echo "-----------------------------------------------------------"
echo ""
echo ""
