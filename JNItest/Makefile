program_NAME 					:= pairhmm
program_C_SRCS 				:= $(wildcard *.c)
program_C_OBJS 				:= ${program_C_SRCS:.c=.o}
program_OBJS 					:= $(program_C_OBJS)
program_INCLUDE_DIRS 	:= ./libcxl
program_LIBRARY_DIRS 	:= ./libcxl
program_LIBRARIES 		:= 
program_STATIC_LIBS		:= ./libcxl/libcxl.a

CFLAGS 	+= $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir)) 
CFLAGS	+= -std=c99 -Wall -m64 -lrt -lpthread -fopenmp -I/usr/lib/jvm/java-8-oracle/include -I/usr/lib/jvm/java-8-oracle/include/linux -fPIC  -MMD -MP -MF "/home/bampetas/Documents/JNI test/pairhmm.o.d"
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

.PHONY: all clean distclean

all: $(program_NAME)

debug: CFLAGS += -DDEBUG
debug: $(program_OBJS)
	gcc  $(program_NAME) $(program_OBJS) $(program_STATIC_LIBS) $(LDFLAGS) $(CFLAGS)
	
$(program_NAME): $(program_OBJS)
	gcc  $(program_NAME) $(program_OBJS) $(program_STATIC_LIBS) $(LDFLAGS) $(CFLAGS)

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)

distclean: clean

