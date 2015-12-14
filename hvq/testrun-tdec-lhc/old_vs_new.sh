rm -r plots
mkdir plots
cp newgenplots.sh plots


# (1) Run old Sudakov
./clean.sh;
cd ..
make clean
make -j16
make -j16 lhef_analysis
#make -j16 main-PYTHIA8-lhef
cd testrun-tdec-lhc
../pwhg_main
cp pwgevents.lhe events/pwgevents_old.lhe
../lhef_analysis
cp pwgLHEF_analysis.top plots/old.top
# ../main-PYTHIA8-lhef
# cp pwgPOWHEG+PYTHIA8-output.top plots/on_PYTHIA_real.top


# (2) get ready to make new sudakov with new seed
#sed -i '' 's/iseed\ 352345/iseed\ 350000/g' powheg.input
sed -i '' 's/newsuda\ 0/newsuda\ 1/g' powheg.input
./clean.sh
cd ..
make clean
make -j16
make -j16 lhef_analysis
#make -j16 main-PYTHIA8-lhef
cd testrun-tdec-lhc
../pwhg_main
cp pwgevents.lhe events/pwgevents_new.lhe
../lhef_analysis
cp pwgLHEF_analysis.top plots/new.top
# ../main-PYTHIA8-lhef
# cp pwgPOWHEG+PYTHIA8-output.top plots/nn_PYTHIA_real.top


# (3) put back to old seed and old sudakov
#sed -i '' 's/iseed\ 350000/iseed\ 352345/g' powheg.input
sed -i '' 's/newsuda\ 1/newsuda\ 0/g' powheg.input


# (4) Plot results
cd plots
./newgenplots.sh *.top
gnuplot *.gp
open *.pdf