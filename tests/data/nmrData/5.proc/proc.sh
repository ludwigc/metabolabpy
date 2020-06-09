#!/bin/tcsh
setenv FID ../50
setenv fidSP fidSP.com
setenv REC2FT recFT.com
setenv in_file nls.in
setenv selection_file nuslist
setenv FST_PNT_PPM 11.5
setenv ROISW 12
setenv proc_out test.dat
setenv SPARSE               y
setenv NUS_TABLE_ORDER      '1'
setenv NUS_TABLE_OFFSET     0
setenv NUS_POINTS           512
setenv NI                   '512 1 '
setenv NIMAX                '2048 1 '
setenv NDIM                 2
setenv MDDTHREADS           16
setenv METHOD               CS
setenv CS_niter             20
setenv CS_alg               IRLS
setenv CS_norm              0
setenv CS_lambda            1.0
setenv CS_VE                n
setenv SRSIZE               0.1

mddnmr4pipeN.sh  1 2 3 4 5
set ecode=$? ;if ( $ecode ) exit($ecode)

