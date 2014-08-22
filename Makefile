parser: parser.cpp
	g++ --std=c++0x parser.cpp -o parser -lnetpbm

clean:
	rm parser
