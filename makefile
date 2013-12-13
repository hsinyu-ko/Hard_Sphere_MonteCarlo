#  FC=gfortran
FC=ifort
FF=-fast
SRC=src/module.f90\
		src/tools.f90\
		src/mc_hard_sphere_nvt.f90
PROGRAM=hard_sp_nvt.x
INP_L=inp.low_density
INP_H=inp.high_density
all: rl rh
$(PROGRAM): clean
	$(FC) $(FF) -o $(PROGRAM) $(SRC)
	rm *.mod
clean:
	rm -f *.o $(PROGRAM) *.mod g_2* r_av restart.dat
veryclean: clean
	rm -f out/* *.png
rl: $(PROGRAM)
	if [ -f restart_l.dat ]; then cp restart_l.dat restart.dat; fi
	./$(PROGRAM) < inp/$(INP_L)
	rm -f ref/g2_ref.dat
	ln -s ref/g2_ref_l.dat ref/g2_ref.dat
	./scripts/plot_rav.py
	./scripts/plot_g2.py
	mv out/r_av{,_l}
	mv out/g_2.dat g2_l
	mv out/restart{,_l}.dat
	mv eq.png eq_l.png
	mv g2.png g2_l.png
rh: $(PROGRAM)
	if [ -f restart_h.dat ]; then cp restart_h.dat restart.dat; fi
	./$(PROGRAM) < inp/$(INP_H)
	rm -f ref/g2_ref.dat
	ln -s ref/g2_ref_h.dat ref/g2_ref.dat
	./scripts/plot_rav.py
	./scripts/plot_g2.py
	mv out/r_av{,_h}
	mv out/g_2.dat g2_h
	mv out/restart{,_h}.dat
	mv eq.png eq_h.png
	mv g2.png g2_h.png
