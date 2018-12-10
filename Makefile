CXX=g++ -g -std=c++0x
LDFLAGS = -I. -lm -O2 -lpthread
PREFIX    = $(DESTDIR)/usr/local
#PREFIX = /usr/local

OBJDIR=.

SRCS = $(wildcard filtering/*.cc)
OBJS = $(foreach OBJ,$(SRCS:.cc=.o),$(OBJDIR)/$(OBJ))
	DEPS = $(wildcard *.h)

$(OBJDIR):
		if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi

$(OBJDIR)/%.o: %.cc $(DEPS)
		$(CXX) -c $(LDFLAGS) -o $@ $<

all: $(OBJDIR) adVNTR-Filtering

adVNTR-Filtering: $(OBJS)
		$(CXX) -o $@ $^ $(LDFLAGS)

.PHONY: clean
.PHONY: all
.PHONY: archive
.PHONY: install
.PHONY: uninstall

clean:
		rm -f *~ $(OBJDIR)/*.o filtering/*.o adVNTR-Filtering

archive: clean

install: adVNTR-Filtering
		install -m 755 adVNTR-Filtering $(DESTDIR)$(PREFIX)/bin

uninstall:
		rm -f $(DESTDIR)$(PREFIX)/bin/adVNTR-Filtering

