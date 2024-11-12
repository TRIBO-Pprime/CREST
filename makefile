include makefile.inc

MODDIR = mod
OBJDIR = obj
SRCDIR = src
INCDIR = inc
LIBDIR = lib

ALL_LIBS  = $(addprefix $(fftw3_dir)/$(LIBDIR)/, *.a)

SGALIB = -Wl,--start-group $(ALL_LIBS) -Wl,--end-group

EXT = f90

VPATH = $(MODDIR):$(OBJDIR):$(SRCDIR)

EXE = main

SRC = $(notdir $(wildcard $(SRCDIR)/*.$(EXT)))

OBJ       = $(SRC:.$(EXT)=.o)
ALL_OBJS  = $(addprefix $(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(algen_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(chole_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(digis_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(fftw3_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(files_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(gplot_dir)/$(OBJDIR)/, *.o)
#~ ALL_OBJS += $(addprefix $(intpl_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(least_dir)/$(OBJDIR)/, *.o)
#~ ALL_OBJS += $(addprefix $(minim_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(qsort_dir)/$(OBJDIR)/, *.o)
#~ ALL_OBJS += $(addprefix $(splin_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(tchev_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(utils_dir)/$(OBJDIR)/, *.o)

ALL_OBJS += $(addprefix $(filtr_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(stats_dir)/$(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(aniso_dir)/$(OBJDIR)/, *.o)

ALL_MODS  = -I$(MODDIR)
ALL_MODS += -I$(algen_dir)/$(MODDIR)
ALL_MODS += -I$(chole_dir)/$(MODDIR)
ALL_MODS += -I$(digis_dir)/$(MODDIR)
ALL_MODS += -I$(fftw3_dir)/$(MODDIR)
ALL_MODS += -I$(files_dir)/$(MODDIR)
ALL_MODS += -I$(gplot_dir)/$(MODDIR)
#~ ALL_MODS += -I$(intpl_dir)/$(MODDIR)
ALL_MODS += -I$(least_dir)/$(MODDIR)
#~ ALL_MODS += -I$(minim_dir)/$(MODDIR)
ALL_MODS += -I$(qsort_dir)/$(MODDIR)
#~ ALL_MODS += -I$(splin_dir)/$(MODDIR)
ALL_MODS += -I$(tchev_dir)/$(MODDIR)
ALL_MODS += -I$(utils_dir)/$(MODDIR)

ALL_MODS += -I$(filtr_dir)/$(MODDIR)
ALL_MODS += -I$(stats_dir)/$(MODDIR)
ALL_MODS += -I$(aniso_dir)/$(MODDIR)

CFLAGS  = $(ALL_MODS) -fopenmp
CFLAGS += -cpp -ffree-form -ffree-line-length-none -march=native -fimplicit-none -fall-intrinsics -fmax-errors=10 -finit-real=nan -ffpe-summary=none

LFLAGS  = $(SGALIB)
LFLAGS += -lpthread -lm -lgomp

ifeq ($(DEBUG),TRUE)
	CFLAGS += -Og -g -Wall -Wextra -fbacktrace -pedantic -fbounds-check -Wuninitialized -fimplicit-none
else
	CFLAGS += -$(OPTC)
endif

ifeq ($(GPROF),TRUE)
	CFLAGS += -pg -g
	LFLAGS += -pg
endif

%.o:	%.$(EXT)
	$(FORT) $(CFLAGS) -c $< -o $(OBJDIR)/$@
	@find . -maxdepth 1 -name '*.mod*' -type f -print0 | xargs -0r mv -t ./$(MODDIR)

$(EXE):	$(OBJ)
	$(FORT) $(ALL_OBJS) $(LFLAGS) -o $(EXE)

mod_crest_param.o :
mod_func_acf.o : mod_crest_param.o
mod_skku_profiles.o : mod_crest_param.o
mod_script.o : mod_skku_profiles.o mod_func_acf.o

main.o : mod_script.o

.PHONY: clean debug gprof all

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(MODDIR)/*.mod
	rm -f $(EXE)

debug:
	make "DEBUG=TRUE"

gprof:
	make "GPROF=TRUE" "DEBUG=TRUE"

all:
	(make clean ; make)


