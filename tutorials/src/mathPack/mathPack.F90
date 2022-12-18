!! ---------------------------------------------------------------------
Module modMathPack
!! ---------------------------------------------------------------------

!! ---------------------------------------------------------------------
Implicit None
!! ---------------------------------------------------------------------

!!... Subroutine Interface
Interface InvMatrix
    Procedure :: InvMatrix_RealReal
    Procedure :: InvMatrix_ComplexComplex
End Interface

!! ---------------------------------------------------------------------
Contains
!! ---------------------------------------------------------------------

    Subroutine InvMatrix_RealReal( A, b, x )
        Implicit None
        Real(8), allocatable :: A(:, :), b(:, :)
        Real(8), allocatable :: x(:, :)
        Integer :: nRow, nCol, nRowB, nColB

        Write(*,*) "Subroutine: InvMatrix_RealReal"

        nRow = size(A, 1)
        nCol = size(A, 2)

        nRowB = size(b, 1)
        nColB = size(b, 2)

        !!... Allocate the solution matrix
        If (Allocated(x)) Deallocate(x)

        BLOCK
            Integer :: M, N, NRHS, NB, MN, LWORK
            Real(8), Allocatable :: WORK(:)
            Integer :: INFO
            Integer, External :: ILAENV

            !!... Match Information
            M    = nRow
            N    = nCol
            NRHS = nColB

            !!... Optimum Block Size (from Lapack driver, DGELS)
            NB = ILAENV(1, 'DGEQRF', '', M, N, -1, -1)

            !!... From explanation of Lapack package
            MN = min(M, N)
            LWORK = max( 1, MN + max( MN, NRHS )*NB )
            If (Allocated(Work)) Deallocate(Work)
            Allocate( Work(LWORK) )

            !!... Inverse of Real Matrix (from Lapack)
            Call DGELS("N", nRow, nCol, nColB, &
            &   A, nRow, b, nRowB, WORK, LWORK, INFO )

        End Block

        !!... Allocate the solution
        allocate(x, source = b)

    End Subroutine

    Subroutine InvMatrix_ComplexComplex( A, b, x )
        Implicit None
        Complex(8), dimension(:, :) :: A, b
        Complex(8), dimension(:, :) :: x

        Write(*,*) "Subroutine: InvMatrix_ComplexComplex"

        !! Try it !

    End Subroutine

!! ---------------------------------------------------------------------
End Module
!! ---------------------------------------------------------------------
