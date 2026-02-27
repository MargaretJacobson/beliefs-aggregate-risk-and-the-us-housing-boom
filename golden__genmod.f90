        !COMPILER-GENERATED INTERFACE MODULE: Sun Jan 11 17:10:51 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GOLDEN__genmod
          INTERFACE 
            RECURSIVE FUNCTION GOLDEN(AX,BX,CX,F,TOL,XMIN)
              REAL(KIND=8), INTENT(IN) :: AX
              REAL(KIND=8), INTENT(IN) :: BX
              REAL(KIND=8), INTENT(IN) :: CX
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: XMIN
              REAL(KIND=8) :: GOLDEN
            END FUNCTION GOLDEN
          END INTERFACE 
        END MODULE GOLDEN__genmod
