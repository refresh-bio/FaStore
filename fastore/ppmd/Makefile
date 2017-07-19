## Take a look at PPMdType.h for additional compiler & environment options
LIB = libppmd
DEBUG = 0
CPP_SET = PPMd.cpp Model.cpp
C_SET =

CC = gcc
LINK = gcc
CODE_FLAGS = -m64 -fno-exceptions -fno-rtti -Wall -static \
	-Wno-unknown-pragmas -Wno-sign-compare -Wno-conversion

AR_FLAGS = -cvq

ifeq ($(DEBUG),0)
  OPT_FLAGS = -O1 -funroll-all-loops
  DEBUG_FLAGS = -g0
  LD_FLAGS = -s -O3
else
  OPT_FLAGS = -O0
  DEBUG_FLAGS = -g
  LD_FLAGS = -g
endif

OBJ_SET = $(CPP_SET:.cpp=.o) $(C_SET:.c=.o)

.cpp.o:
	@$(CC) $(CODE_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS) -c $^
.c.o:
	@$(CC) $(CODE_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS) -c $^

lib: $(OBJ_SET)
	ar $(AR_FLAGS) $(LIB).a $(OBJ_SET)

clean:
	rm *.o
	rm $(LIB).a