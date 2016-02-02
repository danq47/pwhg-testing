
# (1) clean and remake everything
#./clean.sh
#cd ..
#make clean
#make -j16
#make -j16 lhef_analysis
#cd testrun-tdec-lhc
rm -r plots
mkdir plots
cp newgenplots.sh plots

# (1.5) Loop old, then new, changing the seed each time

for i in {1..1000}
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
done

# Change seed back to 1 for next run
sed -i "s/iseed\ "$[1+$i]"/iseed\ 1/g" powheg.input



