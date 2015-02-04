all: graphGen.exe

%.exe: %.o
	g++ -o $@ $<

%.o: %.cpp
	g++ -c $< -o $@
