    CC := gfortran
    CFLAGS := -O2
    MODULOS = *.o
    .PHONY : clean
    all : $(MODULOS) %.o : %.c
         $(CC) $(CFLAGS) –c $<.c –o $@
    ventana.o : ventana.c bd.o : bd.c gestion.o : gestion.c ventana.o bd.o
         $(CC) $(CFLAGS) –c $<.c –o $@
         $(CC) $* -o $@
    juego: juego.c ventana.o bd.o gestion.o
         $(CC) $(CFLAGS) –c $<.c –o $@
         $(CC) $* -o $@
    clean:
         rm –f $(MODULOS)