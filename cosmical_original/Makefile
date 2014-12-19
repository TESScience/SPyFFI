CR_IMG3_OBJS = cr_img3.o seed_tw_ran.o twister.o 
COSMICAL_OBJS = cosmical.o seed_tw_ran.o twister.o 

twister.o:	twister.c
		gcc -c -o twister.o twister.c

seed_tw_ran.o:	seed_tw_ran.c
		gcc -c -o seed_tw_ran.o seed_tw_ran.c 

gasdev2t.o:	gasdev2t.c
		gcc -c -o gasdev2t.o gasdev2t.c

cr_img3.o:	cr_img3.c
		gcc -c -o cr_img3.o cr_img3.c

cosmical.o:	cosmical.c
		gcc -c -o cosmical.o cosmical.c


cr_img3:	${CR_IMG3_OBJS}
		gcc -o cr_img3 ${CR_IMG3_OBJS} -lm

cosmical:	${COSMICAL_OBJS}
		gcc -o cosmical ${COSMICAL_OBJS} -lm