
# GCC
CC 	= gcc
CFLAGS1 = -DGCC -Wall
CFLAGS2 = -DGCC -Wall -funroll-loops -funroll-all-loops -O2
CFLAGS3 = -DGCC -Wall -funroll-loops -funroll-all-loops -O3

# Solaris C
# CC 	= /opt/SUNWspro/bin/cc
# CC 	= /export/home/Forte6/SUNWspro/WS6U2/bin/cc
# CFLAGS1 = -g 
# CFLAGS2 = -g -fast
# CFLAGS3 = -g -fast -xarch=v9a


all: sorting sortingO2 sortingO3

sorting: sorting.c
#	${CC} -DDO_ISORT=0 ${CFLAGS1} -o sorting     sorting.c -lm
	${CC} -DDO_ISORT=1 ${CFLAGS1} -o sorting_all sorting.c -lm

sortingO2: sorting.c
#	${CC} -DDO_ISORT=0 ${CFLAGS2} -o sortingO2     sorting.c -lm
	${CC} -DDO_ISORT=1 ${CFLAGS2} -o sortingO2_all sorting.c -lm

sortingO3: sorting.c
#	${CC} -DDO_ISORT=0 ${CFLAGS3} -o sortingO3     sorting.c -lm
	${CC} -DDO_ISORT=1 ${CFLAGS3} -o sortingO3_all sorting.c -lm

clean: 
	rm -f core *~ \
	sorting sortingO2 sortingO3 \
	sorting_all sortingO2_all sortingO3_all
