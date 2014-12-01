C5G7 macroscopic cross section data
 1 3
 2.0E+07 
!
!Comments can appear after the first 3 lines and between macro/micro blocks
!
!In the second line, the first number is number of groups and the other is
!number of cross section sets. 
!
!In the third line the energy group bounds are made up.
! 
!Data here is derived from NEA/NSC/DOC(2003)16 or ISBN 92-64-02139-6
!Table 1 of Appendix A.
!
!The control rod cross sections come from NEA/NSC/DOC(2005)16 or 
!ISBN 92-64-01069-6 Table 1 of Appendix A.
!
!  Abs       nu-fiss       fiss        chi
! scat mat
!UO2 fuel-clad  
XSMACRO UO2-3.3 0
  7.00000E-01 8.00000E-01 3.33333E-01 1.0000E+00
  2.20000E+00


!Guide tube
XSMACRO GuideTube 0
  3.5000E-02 0.000000E+00 0.00000E+00 0.0000E+00
  2.00000E+00


!Moderator
XSMACRO Moderator 0
  6.0000E-02 0.000000E+00 0.00000E+00 0.0000E+00
  4.00000E+00
