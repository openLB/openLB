#  This file is part of the OpenLB library
#
#  Copyright (C) 2007 Mathias Krause
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
## definitions

include definitions.mk

include $(ROOT)/global.mk

OBJECTS    := $(foreach file, $(SRC), $(PWD)/$(file:.cpp=.o))
DEPS       := $(foreach file, $(SRC), $(PWD)/$(file:.cpp=.d))

###########################################################################
## all

all : depend compile link


###########################################################################
## dependencies

depend : $(DEPS)

$(PWD)/%.d : %.cpp
	@echo Create dependencies for $<
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $(IDIR) $< \
                       | sed -e "s!$*\.o!$(PWD)\/$*\.o!1" > .tmpfile; \
                       cp -f .tmpfile $@;'

###########################################################################
## compile

compile : $(OBJECTS)

$(PWD)/%.o: %.cpp
	@echo Compile $<
	$(CXX) $(CXXFLAGS) $(IDIR) -c $< -o $@

###########################################################################
## clean

clean : cleanrub cleanobj cleandep

cleanrub:
	@echo Clean rubbish files
	@rm -f *~ core .tmpfile tmp/*.* $(OUTPUT)
	@rm -f tmp/vtkData/*.* tmp/vtkData/data/*.* tmp/imageData/*.* tmp/gnuplotData/*.* tmp/gnuplotData/data/*.*

cleanobj:
	@echo Clean object files
	@rm -f $(OBJECTS)

cleandep:
	@echo Clean dependencies files
	@rm -f $(DEPS)

cleanbuild:
	@cd $(ROOT); \
	 $(MAKE) cleanlib;

###########################################################################
## update lib

$(ROOT)/$(LIBDIR)/lib$(LIB).a :
	@cd $(ROOT); \
	$(MAKE) all

###########################################################################
## link

link: $(OUTPUT)

$(OUTPUT): $(OBJECTS) $(ROOT)/$(LIBDIR)/lib$(LIB).a
	@echo Link $@
	$(CXX) $(foreach file, $(SRC), $(file:.cpp=.o)) $(LDFLAGS) -L$(ROOT)/$(LIBDIR) -l$(LIB) -lz -o $@

###########################################################################
## include dependencies

ifneq "$(strip $(wildcard *.d))" ""
   include $(foreach file,$(DEPS),$(file))
endif

###########################################################################
###########################################################################
