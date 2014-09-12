
bill: graphics_for_bill.o
	g++ $(CFLAGS) -o bill graphics_for_bill.o -lGL -lm -lgsl -lgslcblas -lglut

graphics_for_bill.o: graphics_for_bill.cpp
	g++ $(CFLAGS) -c graphics_for_bill.cpp

clean:
	rm -f *.o bill







