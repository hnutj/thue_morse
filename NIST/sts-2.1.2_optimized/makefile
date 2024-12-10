CC = /usr/bin/gcc
C++ = /usr/bin/g++
NVCC = /usr/local/cuda-9.0/bin/nvcc
GCCFLAGS = -c -Wall -fopenmp -std=c11
CUDAFLAGS = -L/usr/local/cuda-9.0/lib64 -lcudart -lcufft
ROOTDIR = .
SRCDIR = $(ROOTDIR)/src
OBJDIR = $(ROOTDIR)/obj
VPATH  = src:obj:include

OBJ = $(OBJDIR)/cephes.o \
	  $(OBJDIR)/assess.o $(OBJDIR)/frequency.o $(OBJDIR)/blockFrequency.o \
      $(OBJDIR)/cusum.o $(OBJDIR)/runs.o $(OBJDIR)/longestRunOfOnes.o \
      $(OBJDIR)/serial.o $(OBJDIR)/rank.o $(OBJDIR)/discreteFourierTransform.o \
      $(OBJDIR)/overlappingTemplateMatchings.o $(OBJDIR)/universal.o \
      $(OBJDIR)/approximateEntropy.o $(OBJDIR)/randomExcursions.o \
      $(OBJDIR)/randomExcursionsVariant.o $(OBJDIR)/linearComplexity.o \
      $(OBJDIR)/cufft.o $(OBJDIR)/dfft.o $(OBJDIR)/matrix.o \
      $(OBJDIR)/utilities.o $(OBJDIR)/generators.o $(OBJDIR)/genutils.o \
	  $(OBJDIR)/nonOverlappingHelper.o $(OBJDIR)/interface.o \
	  $(OBJDIR)/nonOverlappingTemplateMatchings.o \
	  $(OBJDIR)/dm.o $(OBJDIR)/jumphop.o \
	  $(OBJDIR)/selfCorrelation.o $(OBJDIR)/binaryDerivate.o \
	  $(OBJDIR)/pokerDetect.o $(OBJDIR)/runsDistribution.o

assess: $(OBJ)
	$(CC) -o $@ $(OBJ) -lm -fopenmp -std=c11 $(CUDAFLAGS)

$(OBJDIR)/assess.o: $(SRCDIR)/assess.c defs.h decls.h utilities.h
	$(CC) -o $@ -c $(SRCDIR)/assess.c

$(OBJDIR)/frequency.o: $(SRCDIR)/frequency.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/frequency.c

$(OBJDIR)/blockFrequency.o: $(SRCDIR)/blockFrequency.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/blockFrequency.c

$(OBJDIR)/cusum.o: $(SRCDIR)/cusum.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/cusum.c

$(OBJDIR)/runs.o: $(SRCDIR)/runs.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/runs.c

$(OBJDIR)/longestRunOfOnes.o: $(SRCDIR)/longestRunOfOnes.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/longestRunOfOnes.c

$(OBJDIR)/rank.o: $(SRCDIR)/rank.c defs.h externs.h matrix.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/rank.c

$(OBJDIR)/cufft.o: \
        $(SRCDIR)/cufft.cu
	$(NVCC) -o $@ -c $(SRCDIR)/cufft.cu -lcufft --compiler-bindir /usr/bin/gcc-4.8.5

$(OBJDIR)/discreteFourierTransform.o: $(SRCDIR)/discreteFourierTransform.c \
        defs.h externs.h utilities.h
	$(CC) -o $@ $(GCCFLAGS) $(FFTFLAGS) $(SRCDIR)/discreteFourierTransform.c

$(OBJDIR)/nonOverlappingHelper.o: \
        $(SRCDIR)/nonOverlappingHelper.cu
	$(NVCC) -o $@ -c $(SRCDIR)/nonOverlappingHelper.cu --compiler-bindir /usr/bin/gcc-4.8.5

$(OBJDIR)/interface.o: $(SRCDIR)/interface.cpp
	$(C++) -o $@ -c $(SRCDIR)/interface.cpp

$(OBJDIR)/nonOverlappingTemplateMatchings.o: \
        $(SRCDIR)/nonOverlappingTemplateMatchings.c defs.h externs.h utilities.h 
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/nonOverlappingTemplateMatchings.c

$(OBJDIR)/overlappingTemplateMatchings.o: \
        $(SRCDIR)/overlappingTemplateMatchings.c defs.h externs.h utilities.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/overlappingTemplateMatchings.c 

$(OBJDIR)/universal.o: $(SRCDIR)/universal.c defs.h externs.h utilities.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/universal.c

$(OBJDIR)/approximateEntropy.o: $(SRCDIR)/approximateEntropy.c defs.h externs.h utilities.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/approximateEntropy.c

$(OBJDIR)/randomExcursions.o: $(SRCDIR)/randomExcursions.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/randomExcursions.c

$(OBJDIR)/randomExcursionsVariant.o: $(SRCDIR)/randomExcursionsVariant.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/randomExcursionsVariant.c

$(OBJDIR)/serial.o: $(SRCDIR)/serial.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/serial.c

$(OBJDIR)/linearComplexity.o: $(SRCDIR)/linearComplexity.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/linearComplexity.c

$(OBJDIR)/dm.o: $(SRCDIR)/dm.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/dm.c

$(OBJDIR)/jumphop.o: $(SRCDIR)/jumphop.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/jumphop.c

$(OBJDIR)/selfCorrelation.o: $(SRCDIR)/selfCorrelation.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/selfCorrelation.c

$(OBJDIR)/binaryDerivate.o: $(SRCDIR)/binaryDerivate.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/binaryDerivate.c

$(OBJDIR)/pokerDetect.o: $(SRCDIR)/pokerDetect.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/pokerDetect.c

$(OBJDIR)/runsDistribution.o: $(SRCDIR)/runsDistribution.c defs.h externs.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/runsDistribution.c

$(OBJDIR)/dfft.o: $(SRCDIR)/dfft.c
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/dfft.c

$(OBJDIR)/matrix.o: $(SRCDIR)/matrix.c defs.h externs.h utilities.h matrix.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/matrix.c

$(OBJDIR)/genutils.o: $(SRCDIR)/genutils.c config.h genutils.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/genutils.c

$(OBJDIR)/cephes.o: $(SRCDIR)/cephes.c cephes.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/cephes.c

$(OBJDIR)/utilities.o: $(SRCDIR)/utilities.c defs.h externs.h utilities.h config.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/utilities.c

$(OBJDIR)/generators.o: $(SRCDIR)/generators.c defs.h externs.h utilities.h \
        config.h generators.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/generators.c

clean:
	rm -f assess $(OBJDIR)/*.o

rebuild: clean assess
