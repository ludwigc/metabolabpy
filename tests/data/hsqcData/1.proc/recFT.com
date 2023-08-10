#!/bin/tcsh -f
echo '| ' in $0 $1
if( $#argv < 1 ) then
echo "Use: $0 <input pipe> <template for output spectrum>"
echo "nmrPipe processing of YZ dimensions after MDD reconstruction"
exit 1
endif
set ft4trec=$1
if( $#argv > 1 ) set proc_out=$2
if( ! -f $ft4trec ) then
ls $ft4trec
echo $0 failed
exit 2
endif
echo '|   Processing time domain MDD reconstruction '
echo
echo Processing Y dimensions
showhdr $ft4trec
cat $ft4trec                                        \
| nmrPipe  -fn TP -auto                             \
| nmrPipe  -fn SP -off 0.48 -end 0.95 -pow 2 -c 0.5 \
| nmrPipe  -fn ZF -size 4096	                    \
| nmrPipe  -fn FT                                   \
| nmrPipe  -fn PS -hdr                              \
| nmrPipe  -fn PS -p0 -6.2 -p1 3.0 -di              \
| nmrPipe  -fn TP  -auto                            \
| nmrPipe  -fn POLY -auto \
| nmrPipe  -fn TP -auto \
| nmrPipe  -fn POLY -auto \
| nmrPipe  -fn TP  -auto                            \
-ov -out $proc_out
echo $proc_out ready
exit
