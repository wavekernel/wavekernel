.PHONY: all clean

all:
	make -C eigenkernel
	cp eigenkernel/src/*.mod eigenkernel/src/libeigenkernel.a ./
	make -C wavepacket
	cp wavepacket/src/*.mod wavepacket/src/libwavepacket.a ./
clean:
	make -C eigenkernel clean
	make -C wavepacket clean
	rm -f ./*.mod ./lib*.a
