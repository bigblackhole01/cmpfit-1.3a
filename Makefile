
OFILES=mpfit.o
LIBFILE=libmpfit.a

all: $(LIBFILE) testmpfit

clean:
	rm -f $(OFILES) testmpfit $(LIBFILE)

mpfit.o: mpfit.c mpfit.h
	clang -c -o $@ $< $(CFLAGS)

$(LIBFILE): $(OFILES)
	$(AR) r $@ $(OFILES)

testmpfit: testmpfit.c libmpfit.a
	clang -o $@ $(CFLAGS) testmpfit.c -L. -lmpfit -O3 -lm -fopenmp 
