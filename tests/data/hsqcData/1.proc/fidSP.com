bruk2pipe -in ./ser -bad 0.0 -noswap -DMX -decim 1600 -dspfvs 21 -grpdly 76 \
-xN       	1024      	-yN       	4096      	 \
-xT       	512       	-yT       	2048      	 \
-xMODE    	Complex   	-yMODE    	Echo-AntiEcho	 \
-xSW      	12500.000 	-ySW      	38167.939 	 \
-xOBS     	799.750   	-yOBS     	201.097   	 \
-xCAR     	4.700     	-yCAR     	80.000    	 \
-xP0      	-78.0     	-yP0      	 90.0     	 \
-xP1      	 -4.0     	-yP1      	  0.0     	 \
-xLAB     	1H        	-yLAB     	13C       	 \
-ndim 2	-aq2D States \
| nmrPipe -fn POLY -time                               \
| nmrPipe -fn SP -off 0.48 -end 0.95 -pow 2 -c 0.5     \
| nmrPipe -fn ZF -auto                                 \
| nmrPipe -fn FT -auto                                 \
| nmrPipe -fn PS -hdr                                  \
| nmrPipe -fn PS -p0 28.0809 -p1 -3.1788 -di                       \
| nmrPipe -fn POLY -auto -xn 5.0ppm -ord 1              \
| nmrPipe -fn EXT -sw                                  \
| pipe2xyz -z -out ft/data%03d.DAT -ov -nofs -verb
