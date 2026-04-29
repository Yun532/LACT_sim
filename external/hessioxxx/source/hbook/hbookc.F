C     +++ SccsId: "@(#)hbookc.F	1.3  12/1/98" +++
C     
C
C@@ Interface functions for using HBOOK from C functions.

C------------------------------------------------------
C@@@BEGINHBOOK()
C@@    Initialize HBOOK interface and set memory limit.
C@
C@     Arguments:  (none)
C@
C@     Return value:  (none)
C@
C------------------------------------------------------

      Subroutine BEGINHBOOK
#ifdef PAWMEMORY
      Common/PAWC/B(PAWMEMORY)
#else
      Common/PAWC/B(20000000)
#endif
      Integer init
      Save init
      Data init /0/
      If ( init.Eq.0 ) Then
         init = 1
#ifdef PAWMEMORY
         Call HLIMIT(PAWMEMORY)
#else
         Call HLIMIT(20000000)
#endif
      Endif
      End

C------------------------------------------------------
C@@@SAVEHBOOK()
C@@    Write all booked HBOOK histograms to a file.
C@
C@     Arguments:  fname  --  HBOOK output file
C@
C@     Return value:  (none)
C@
C------------------------------------------------------

      Subroutine SAVEHBOOK (fname, id)
      Character fname*(*)
      Integer id
      CHARACTER*255 CHSAVE

C... Rather than simply calling
C...       Call HRPUT(0,fname,' ')
C... we use now code from HRPUT to allow case-sensitive
C... file names:

      CALL RZCDIR(CHSAVE,'R')
      LSAVE = LENOCC(CHSAVE)
      CALL HROPEN(88,'HRPUT',fname,'NP',1024,ISTAT)
      IF(ISTAT.NE.0) RETURN
      CALL HROUT(id,ICYCLE,' ')
      CALL HRENDC('HRPUT')
      CLOSE(88)
      CALL RZCDIR(CHSAVE(1:LSAVE),' ')

      End

C------------------------------------------------------
C@@@HBOOKN()
C@@    book a ntuple
C@
C@     Arguments:  fname  --  HBOOK output file
C@
C@     Return value:  (none)
C@
C------------------------------------------------------

      Subroutine INITHBOOKN (id,chtitle,nvar,nprime,tags)
      Character chtitle*(*)
      INTEGER id
      INTEGER nvar
      INTEGER nprime
      Character tags(*)*(*)      

      Call HBOOKN(id,chtitle,nvar,' ',nprime,tags)
      End
