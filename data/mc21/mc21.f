C Copyright (C) 2016 University of Vienna
C All rights reserved.
C BSD license.
C Author: Ali Baharev <ali.baharev@gmail.com>
C
        SUBROUTINE MAXMATCH(N, ICN, LICN, IP, LENR, IPERM, NUMNZ, IW, 
     C                      LIW)
        INTEGER N
        INTEGER ICN(LICN)
        INTEGER LICN
        INTEGER IP(N)
        INTEGER LENR(N)
        INTEGER IPERM(N)
        INTEGER NUMNZ
        INTEGER IW(LIW)
Cf2py integer, intent(in) :: N
Cf2py integer, intent(in) :: ICN(LICN)
Cf2py integer, intent(in) :: LICN
Cf2py integer, intent(in) :: IP(N)
Cf2py integer, intent(in) :: LENR(N)
Cf2py integer, intent(inout) :: IPERM(N)
Cf2py integer, intent(inout) :: NUMNZ
Cf2py integer, intent(inout) :: IW(LIW)
C
C Compile as f2py -c -m mc21 mc21.f
C
        EXTERNAL MC21AD
        CALL MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
        END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 Copy the souce code of MC21AD here! Compile as f2py -c -m mc21 mc21.f