C**************************************************************
C  For calculating bond length
C**************************************************************
      SUBROUTINE BONDS (NATOM,XX,RR)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RR(36,36),XX(36,3)
      
      DO I=1,NATOM
        DO J=1,I
          RR(J,I)=0
        END DO
      END DO

      DO K=1,3
        DO I=1,NATOM
          DO J=1,I
            RR(J,I)=RR(J,I)+(XX(I,K)-XX(J,K))*(XX(I,K)-XX(J,K))
          END DO
        END DO
      END DO
      DO I=1,NATOM
        DO J=1,I
          RR(J,I)=SQRT(RR(J,I))
          RR(I,J)=RR(J,I)
        END DO
      END DO

9200  FORMAT (f8.0,7f10.5)

      RETURN
      END
C**************************************************************
C**************************************************************
