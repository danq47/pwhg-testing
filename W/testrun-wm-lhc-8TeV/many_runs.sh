
# (1) clean and remake everything
./clean.sh
#cd ..
#make clean
#make -j16
#make -j16 lhef_analysis
#cd testrun-tdec-lhc
./clean.sh
rm -r plots
mkdir plots
cp newgenplots.sh plots
cp merge plots
mkdir plots/merged

# (1.5) Loop old, then new, changing the seed each time

# set seed back to 1

sed -i "s/iseed\ .*/iseed\ 1/g" powheg.input
sed -i "s/newsuda\ 1/newsuda\ 0/g" powheg.input


for i in {1..10000}
do
# (2) Run old Sudakov
#	./clean.sh
        rm pwgevents.lhe
	../pwhg_main
	../lhef_analysis
	cp pwgLHEF* plots/old"$i".top


# (3) Change to new Sudakov
	sed -i  "s/newsuda\ 0/newsuda\ 1/g" powheg.input
#	./clean.sh
	rm pwgevents.lhe
	../pwhg_main
	../lhef_analysis
	cp pwgLHEF* plots/new"$i".top

# change back to old suda	
	sed -i "s/newsuda\ 1/newsuda\ 0/g" powheg.input

# change seed for next run (here they simply run with seed=1,seed=2...)
# Would be better to somehow get it to read the seed and then increment it every time
# This would mean that it wouldn't need to start from 1 each time (less room for error)
	sed -i "s/iseed\ "$i"/iseed\ "$[1+$i]"/g" powheg.input

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

# Change seed back to 1 for next run
sed -i "s/iseed\ "$[1+$i]"/iseed\ 1/g" powheg.input



