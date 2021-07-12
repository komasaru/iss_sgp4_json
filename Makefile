gcc_options = -std=c++17 -Wall -O2 --pedantic-errors

iss_sgp4_json: iss_sgp4_json.o eop.o sgp4.o tle.o blh.o time.o
	g++ $(gcc_options) -o $@ $^

iss_sgp4_json.o : iss_sgp4_json.cpp
	g++ $(gcc_options) -c $<

eop.o : eop.cpp
	g++ $(gcc_options) -c $<

sgp4.o : sgp4.cpp
	g++ $(gcc_options) -c $<

tle.o : tle.cpp
	g++ $(gcc_options) -c $<

blh.o : blh.cpp
	g++ $(gcc_options) -c $<

time.o : time.cpp
	g++ $(gcc_options) -c $<

run : iss_sgp4_json
	./iss_sgp4_json

clean :
	rm -f ./iss_sgp4_json
	rm -f ./*.o

.PHONY : run clean

