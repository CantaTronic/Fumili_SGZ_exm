
FDDIR := $(FDMODULE)

CXXFLAGS := -O2 -Wall -fPIC -g
CXXFLAGS += -I$(FDDIR)/include
CXXFLAGS += `root-config --cflags`

LDFLAGS := `root-config --libs` -lMinuit -lGeom
LDFLAGS += -g -O2
LDFLAGS += -L$(FDDIR)/lib -lGEvent -lFDTracks

OBJS := Vika VikaC

all: $(OBJS)

%: %.o
	@echo 'Linking executable $@'
	@$(CXX) $^ -o $@ $(LDFLAGS)

%.o: %.cc %.d
	@echo 'Compiling $@'
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(CPPFLAGS)

clean:
	$(RM) $(OBJS) $(addsuffix .o,$(OBJS)) $(addsuffix .d,$(OBJS))

%.d: %.cc
	@echo Making dependency for file $< ...
	@(echo -n $@ $(dir $@); $(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $<) >>$@

.PHONY: all clean

.DELETE_ON_ERROR:

-include $(addsuffix .d,$(OBJS))
