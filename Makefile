#  This file is part of the OpenLB library
#
#  Copyright (C) 2017 Markus Mohrhard
#  E-mail contact: info@openlb.net
#  The most recent release of OpenLB can be downloaded at
#  <http://www.openlb.net/>
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#  Boston, MA  02110-1301, USA.

###########################################################################

## installing the configuration makefile

ROOT := .

include global.mk

CLEANTARGETS :=

SAMPLESCLEAN :=

LIB_OBJECTS := 

SAMPLES :=

TESTS :=

BENCHMARKS :=

EXTRA_IDIR :=

DEPS := $(shell find $(DEPENDDIR)/*/ -name *.d)

INCLUDEDIR := $(foreach d,$(INCLUDEDIRS),-I./$(d))

###########################################################################
## Building libolb.a and the  

all: compile

include $(addsuffix /module.mk, $(SUBDIRS))

include $(addsuffix /module.mk, $(EXAMPLEDIRS))

include $(addsuffix /module.mk, $(TESTDIRS))

include $(foreach dependency_file,$(LIB_OBJECTS:.o=.d), $(wildcard $(subst $(OBJDIR)/,$(DEPENDDIR)/,$(dependency_file))))

compile: lib

$(OBJDIR)/%.o: %.cpp | $(DEPENDDIR)/%.d
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -DTIXML_USE_STL $(INCLUDEDIR) -c $< -o $@

lib: $(LIBDIR)/lib$(LIB).a $(LIBDIR)/libz.a

$(LIBDIR)/lib$(LIB).a: $(LIB_OBJECTS)
	@mkdir -p $(LIBDIR)
	@$(ARPRG) -rusv $(LIBDIR)/lib$(LIB).a $(LIB_OBJECTS) > /dev/null

.PRECIOUS: $(DEPENDDIR)/%.d

$(DEPENDDIR)/%.d : %.cpp
	@mkdir -p $(dir $@)
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $(INCLUDEDIR) $(EXTRA_IDIR) $< \
	    | sed -e "s!$(basename $(notdir $@))\.o!$(OBJDIR)/$*.o!1" > $@;'

###########################################################################
## removing generated files

cleansamples: $(SAMPLESCLEAN)

cleanbuild cleanlib:
	@rm -f $(LIB_OBJECTS) &> /dev/null || true
	@rm -f $(LIBDIR)/lib$(LIB).a &> /dev/null || true

clean: $(CLEANTARGETS) cleanlib
	@rm -r -f $(shell find $(DEPENDDIR) -name *.d) &> /dev/null || true


## handling of the examples

samples: lib $(SAMPLES)

############################yy###############################################
## user guide documentation

userguide:
	@cd doc/userGuide/; \
	latexmk -pdf -silent -f olb-ug.tex

###########################################################################
## doxygen documentation

doxygen:
	doxygen doc/DoxygenConfig

###########################################################################
## checking whether to automatically install the git hook
GITDEPEND := $(shell test -d ".git" && echo .git/hooks/pre-commit)

ifneq ($(strip $(GITDEPEND)),)
.git/hooks/pre-commit: hooks/pre-commit
	@echo Installing git hook
	cp hooks/pre-commit .git/hooks/pre-commit
	chmod +x .git/hooks/pre-commit
endif
