      PROGRAM ASTECDIAMOND

            DIMENSION IAD(-200:200,-200:200), IVAD(-400:400,-400:400)
            INTEGER dim
            COMMON /AD/ IAD, IVAD, dim

            OPEN(UNIT=1, FILE='diamond.txt', STATUS='UNKNOWN')
            OPEN(UNIT=2, FILE='diamond2.txt', STATUS='UNKNOWN')

            dim = 0

            CAll InitDiamond()

            DO I = 1, 150
                  CALL UpdateDiamond()
            END DO

c           in active cells
            DO I = -dim+1, dim, 2
            DO J = -dim+1, dim, 2

                  IF (IAD(I,J) .EQ. 0) THEN
                        CALL InitSquare(I, J)
                  END IF

            END DO
            END DO

            DO I = -dim+1, dim+1
                  WRITE(1,*) IVAD(I,-dim-1:dim+1)
            END DO

      END PROGRAM ASTECDIAMOND

      SUBROUTINE InitDiamond()
            DIMENSION IAD(-200:200,-200:200), IVAD(-400:400,-400:400)
            COMMON /Ad/ IAD, IVAD, dim

c           init all cells with 0
            DO I = -100, 100
                  DO J = -100, 100
                        IAD(I,J) = 0
                  END DO
            END DO

            DO I = -200, 200
                  DO J = -200, 200
                        IVAD(I,J) = 0
                  END DO
            END DO


      END SUBROUTINE InitDiamond

      SUBROUTINE UpdateDiamond()
            DIMENSION IAD(-200:200,-200:200),IaVAD(-200:200,-200:200),
     +      IVAD(-400:400,-400:400)
            INTEGER dim
            COMMON /AD/ IAD, IVAD, dim

c           in active cells
            DO I = -dim+1, dim, 2
            DO J = -dim+1, dim, 2

                  IF (IAD(I,J) .EQ. 0) THEN
                        CALL InitSquare(I, J)
                  END IF

            END DO
            END DO

c           Update grid size
            dim = dim + 1

c           in active cells
            DO I = -dim+1, dim, 2
            DO J = -dim+1, dim, 2
                  
                  IF (IAD(I,J) .GT. 1) THEN
                        CALL EraseSquare(I, J)
                  END IF

            END DO
            END DO

c           equalize aux grid to zero
            DO I = -dim-10, dim+10
            DO J = -dim-10, dim+10
                  IaVAD(I,J) = 0
            END DO
            END DO

c           Move vertices
            DO I = -dim+1, dim
            DO J = -dim+1, dim

                  IF (IVAD(I,J) .EQ. 1) THEN

                        IaVAD(I+1,J+1) = 1

                        IAD(I+1,J+1) = IAD(I+1,J+1) + 1
                        IAD(I-1,J-1) = IAD(I-1,J-1) - 1

                  ELSE IF (IVAD(I,J) .EQ. 2) THEN

                        IaVAD(I-1,J+1) = 2

                        IAD(I-2,J+1) = IAD(I-2,J+1) + 1
                        IAD(I,J-1) = IAD(I,J-1) - 1

                  ELSE IF (IVAD(I,J) .EQ. -1) THEN

                        IaVAD(I-1,J-1) = -1

                        IAD(I-2,J-2) = IAD(I-2,J-2) + 1
                        IAD(I,J) = IAD(I,J) - 1
                        
                  ELSE IF (IVAD(I,J) .EQ. -2) THEN

                        IaVAD(I+1,J-1) = -2

                        IAD(I+1,J-2) = IAD(I+1,J-2) + 1
                        IAD(I-1,J) = IAD(I-1,J) - 1

                  END IF
            END DO
            END DO

c           copy aux grid to main grid
            DO I = -dim-10, dim+10
            DO J = -dim-10, dim+10
                  IVAD(I,J) = IaVAD(I,J)
            END DO
            END DO

      END SUBROUTINE UpdateDiamond

      SUBROUTINE InitSquare(i, j)
            DIMENSION IAD(-200:200,-200:200), IVAD(-400:400,-400:400)
            COMMON /AD/ IAD, IVAD, dim


            toss = RAND(0)

            IF (toss .LT. 0.5) THEN
                  IVAD(i,j) = -1
                  IVAD(i+1,j+1) = 1

                  IAD(i-1,j-1) = IAD(i-1,j-1) + 1
                  IAD(i+1,j+1) = IAD(i+1,j+1) + 1
            ELSE
                  IVAD(i+1,j) = -2
                  IVAD(i,j+1) = 2

                  IAD(i-1,j+1) = IAD(i-1,j+1) + 1
                  IAD(i+1,j-1) = IAD(i+1,j-1) + 1
            END IF

            IAD(i,j) = IAD(i,j) + 2

      END SUBROUTINE InitSquare

      SUBROUTINE EraseSquare(i, j)
            DIMENSION IAD(-200:200,-200:200), IVAD(-400:400,-400:400)
            COMMON /AD/ IAD, IVAD, dim

            DO k = i, i+1
            DO l = j, j+1
                  IF (IVAD(k,l) .NE. 0) THEN
                        
                  IF ((IVAD(k,l) .EQ. 1) .OR. (IVAD(k,l) .EQ. -1)) THEN
                        IAD(k-1,l-1) = IAD(k-1,l-1) - 1
                        IAD(k,l) = IAD(k,l) - 1
                  ELSE IF ((IVAD(k,l) .EQ. 2).OR.(IVAD(k,l).EQ.-2)) THEN
                        IAD(k-1,l) = IAD(k-1,l) - 1
                        IAD(k,l-1) = IAD(k,l-1) - 1
                  END IF

                  IVAD(k,l) = 0
                  
                  END IF
            END DO
            END DO

      END SUBROUTINE EraseSquare
