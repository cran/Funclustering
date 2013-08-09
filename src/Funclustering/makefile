include $(R_HOME)/etc${R_ARCH}/Makeconf
MYCXXFLAGS = -DNDEBUG -DRPACKAGE
LIB = libFunclustering.a
SRCS = $(wildcard src/IAlgo/*.cpp\
                  src/IModel/*.cpp\
                  src/IStrategy/*.cpp)

OBJS= $(SRCS:%.cpp=%.o)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) ${CPICFLAGS} $(MYCXXFLAGS) $< -c -o $@
	
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -rc $@ $?
  
.PHONY: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o
