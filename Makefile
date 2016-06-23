TARGET = libFitTools.so
SRCFILES = src/TFractionFitterSimul.cc src/TemplateChi2Fitter.cc src/IntegralChi2Fitter.cc src/EfficiencyFitter.cc
HDRFILES = $(patsubst src/%.cc,include/%.h,$(SRCFILES))
ALLSOURCES = $(SRCFILES) src/Dict.cc
ALLOBJECTS = $(patsubst src/%.cc,obj/%.o,$(ALLSOURCES))

CFLAGS = -c -O3 -Wall -fPIC $(shell root-config --cflags)
LFLAGS = -shared

EXTRAINC = -I. -I$(ROOFITSYS)/include
LIBS = $(shell root-config --libs)

all: $(TARGET)

clean:
	rm -rf $(TARGET) obj src/Dict.* Dict_rdict.pcm > /dev/null 2>&1

$(TARGET): $(ALLOBJECTS)
	g++ $(LFLAGS) -o $(TARGET) $(LIBS) $^

obj/Dict.o: src/Dict.cc
	-@mkdir -p obj 2>/dev/null
	g++ $(CFLAGS) $(EXTRAINC) -o $@ $< $(LIBS)

obj/%.o: src/%.cc include/%.h
	-@mkdir -p obj 2>/dev/null
	g++ $(CFLAGS) $(EXTRAINC) -o $@ $< $(LIBS)

src/Dict.cc: $(HDRFILES) include/LinkDef.h
	rootcling -f $@ $^
	mv src/Dict_rdict.pcm .
