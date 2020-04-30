#................................................
#      How to compile whole program 
#................................................

all:$
        cd src; make all; cp *.x ../exe/.
clean: $
        rm exe/*.x; cd src; make clean

