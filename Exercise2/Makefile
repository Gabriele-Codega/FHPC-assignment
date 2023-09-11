### MKL libraries
###
###
MKL= -L${MKLROOT}/lib/intel64  -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

### OpenBLAS libraries 
OPENBLASROOT=${OPENBLAS_ROOT}
### BLIS library
BLISROOT=/u/dssc/gcodeg00/myblis-intel

# precision flag
PFLAG = -DUSE_DOUBLE

cpu: gemm_mkl.x gemm_oblas.x gemm_blis.x


gemm_mkl.x: gemm.c
	gcc $(PFLAG) -DMKL $^ -m64 -I${MKLROOT}/include $(MKL)  -o $@

gemm_oblas.x: gemm.c
	gcc $(PFLAG) -DOPENBLAS $^ -m64 -I${OPENBLASROOT}/include -L/${OPENBLASROOT}/lib -lopenblas -lpthread -o $@ -fopenmp

gemm_blis.x: gemm.c
	gcc $(PFLAG)  -DBLIS $^ -m64 -I${BLISROOT}/include/blis -L/${BLISROOT}/lib -o $@ -lpthread  -lblis -fopenmp -lm
clean:
	rm -rf *.x
