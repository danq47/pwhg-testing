
# (1) clean everything
./clean.sh
rm -r plots
rm merge
mkdir plots
cp newgenplots.sh plots

gfortran -o merge mergedata.f
cp merge plots

mkdir plots/merged

# (1.5) Loop old, then new, changing the seed each time

# set seed back to 1

sed -i "s/iseed\ .*/iseed\ 1/g" powheg.input
sed -i "s/newsuda\ 1/newsuda\ 0/g" powheg.input


for i in {1..10000}
do
# (2) Run old Sudakov
    rm pwgevents.lhe
	../pwhg_main
	../lhef_analysis
	cp pwgLHEF* plots/old"$i".top


# (3) Change to new Sudakov
	sed -i  "s/newsuda\ 0/newsuda\ 1/g" powheg.input
	rm pwgevents.lhe
	../pwhg_main
	../lhef_analysis
	cp pwgLHEF* plots/new"$i".top

# change back to old suda	
	sed -i "s/newsuda\ 1/newsuda\ 0/g" powheg.input

# change seed for next run (here they simply run with seed=1,seed=2...)
	sed -i "s/iseed\ "$i"/iseed\ "$[1+$i]"/g" powheg.input

# Every 25 runs (of both old and new) merge the analysis files to
# stop them from filling up disk space
	if [ $(($i % 25)) -eq 0 ]; then
	    cd plots

	    ./merge 1 old*top
	    mv fort.12 merged/old"$i"-merged.top

	    ./merge 1 new*top
	    mv fort.12 merged/new"$i"-merged.top

	    rm *.top
	    cd ..
	fi

done




