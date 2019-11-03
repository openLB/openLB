# Makefile for zlib
# Copyright (C) 1995-2017 Jean-loup Gailly, Mark Adler
# For conditions of distribution and use, see copyright notice in zlib.h

current_dir := $(dir $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST)))

include global.mk

ZLIB_CFLAGS=$(filter-out -std=%,$(CXXFLAGS))

AR:=$(ARPRG)
ARFLAGS:=rc

ZLIB_SRCDIR:=$(current_dir)
ZLIB_OBJDIR:=$(OBJDIR)/src/external/zlib/
ZINC=

OBJZ := adler32.o crc32.o deflate.o infback.o inffast.o inflate.o inftrees.o trees.o zutil.o
OBJG := compress.o uncompr.o gzclose.o gzlib.o gzread.o gzwrite.o
OBJC := $(foreach file,$(OBJZ) $(OBJG),$(addprefix $(ZLIB_OBJDIR),$(file)))

$(ZLIB_OBJDIR):
	mkdir -p $(ZLIB_OBJDIR)

$(LIBDIR)/libz.a: $(OBJC)
	$(AR) $(ARFLAGS) $@ $(OBJC)

$(ZLIB_OBJDIR)adler32.o: $(ZLIB_SRCDIR)adler32.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)adler32.c

$(ZLIB_OBJDIR)crc32.o: $(ZLIB_SRCDIR)crc32.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)crc32.c

$(ZLIB_OBJDIR)deflate.o: $(ZLIB_SRCDIR)deflate.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)deflate.c

$(ZLIB_OBJDIR)infback.o: $(ZLIB_SRCDIR)infback.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)infback.c

$(ZLIB_OBJDIR)inffast.o: $(ZLIB_SRCDIR)inffast.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)inffast.c

$(ZLIB_OBJDIR)inflate.o: $(ZLIB_SRCDIR)inflate.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)inflate.c

$(ZLIB_OBJDIR)inftrees.o: $(ZLIB_SRCDIR)inftrees.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)inftrees.c

$(ZLIB_OBJDIR)trees.o: $(ZLIB_SRCDIR)trees.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)trees.c

$(ZLIB_OBJDIR)zutil.o: $(ZLIB_SRCDIR)zutil.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)zutil.c

$(ZLIB_OBJDIR)compress.o: $(ZLIB_SRCDIR)compress.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)compress.c

$(ZLIB_OBJDIR)uncompr.o: $(ZLIB_SRCDIR)uncompr.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)uncompr.c

$(ZLIB_OBJDIR)gzclose.o: $(ZLIB_SRCDIR)gzclose.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)gzclose.c

$(ZLIB_OBJDIR)gzlib.o: $(ZLIB_SRCDIR)gzlib.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)gzlib.c

$(ZLIB_OBJDIR)gzread.o: $(ZLIB_SRCDIR)gzread.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)gzread.c

$(ZLIB_OBJDIR)gzwrite.o: $(ZLIB_SRCDIR)gzwrite.c | $(ZLIB_OBJDIR)
	$(CC) $(ZLIB_CFLAGS) $(ZINC) -c -o $@ $(ZLIB_SRCDIR)gzwrite.c


zlib-clean:
	@rm $(ZLIB_OBJDIR)/* &> /dev/null || true
	@rm $(LIBDIR)/libz.a &> /dev/null || true

CLEANTARGETS += zlib-clean

$(ZLIB_OBJDIR)adler32.o $(ZLIB_OBJDIR)zutil.o: $(ZLIB_SRCDIR)zutil.h $(ZLIB_SRCDIR)zlib.h $(ZLIB_SRCDIR)zconf.h
$(ZLIB_OBJDIR)gzclose.o $(ZLIB_OBJDIR)gzlib.o $(ZLIB_OBJDIR)gzread.o $(ZLIB_OBJDIR)gzwrite.o: $(ZLIB_SRCDIR)zlib.h $(ZLIB_SRCDIR)zconf.h $(ZLIB_SRCDIR)gzguts.h
$(ZLIB_OBJDIR)compress.o $(ZLIB_OBJDIR)example.o $(ZLIB_OBJDIR)minigzip.o $(ZLIB_OBJDIR)uncompr.o: $(ZLIB_SRCDIR)zlib.h $(ZLIB_SRCDIR)zconf.h
$(ZLIB_OBJDIR)crc32.o: $(ZLIB_SRCDIR)zutil.h $(ZLIB_SRCDIR)zlib.h $(ZLIB_SRCDIR)zconf.h $(ZLIB_SRCDIR)crc32.h
$(ZLIB_OBJDIR)deflate.o: $(ZLIB_SRCDIR)deflate.h $(ZLIB_SRCDIR)zutil.h $(ZLIB_SRCDIR)zlib.h $(ZLIB_SRCDIR)zconf.h
$(ZLIB_OBJDIR)infback.o $(ZLIB_OBJDIR)inflate.o: $(ZLIB_SRCDIR)zutil.h $(ZLIB_SRCDIR)zlib.h $(ZLIB_SRCDIR)zconf.h $(ZLIB_SRCDIR)inftrees.h $(ZLIB_SRCDIR)inflate.h $(ZLIB_SRCDIR)inffast.h $(ZLIB_SRCDIR)inffixed.h
$(ZLIB_OBJDIR)inffast.o: $(ZLIB_SRCDIR)zutil.h $(ZLIB_SRCDIR)zlib.h $(ZLIB_SRCDIR)zconf.h $(ZLIB_SRCDIR)inftrees.h $(ZLIB_SRCDIR)inflate.h $(ZLIB_SRCDIR)inffast.h
$(ZLIB_OBJDIR)inftrees.o: $(ZLIB_SRCDIR)zutil.h $(ZLIB_SRCDIR)zlib.h $(ZLIB_SRCDIR)zconf.h $(ZLIB_SRCDIR)inftrees.h
$(ZLIB_OBJDIR)trees.o: $(ZLIB_SRCDIR)deflate.h $(ZLIB_SRCDIR)zutil.h $(ZLIB_SRCDIR)zlib.h $(ZLIB_SRCDIR)zconf.h $(ZLIB_SRCDIR)trees.h
