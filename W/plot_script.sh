rm -r nlo_plots
mkdir nlo_plots

cd nlo_plots
cp ../testrun-wm-lhc-8TeV/newgenplots.sh .
cp ../testrun-wm-lhc-8TeV/merge .

cp ../testrun-wm-lhc-8TeV/plots/run* .
cp ../test2/plots/run* .
cp ../test3/plots/run* .

./merge 1 run*top
rm run*top
mv fort.12 LHEF.top

cp ../testrun-wm-lhc-8TeV/pwg-NLO* NLO1.top
cp ../test2/pwg-NLO* NLO2.top
cp ../test3/pwg-NLO* NLO3.top

./merge 1 NLO*top
rm NLO*top
mv fort.12 NLO.top

./newgenplots.sh NLO.top LHEF.top
gnuplot *.gp
open *.pdf
