!> \file m_empire_api.f90

!-----------------------------------------------------------------------
!> \brief The Fortran EMPIRE API interface module.
!> This module must be included in the compile process of a Fortran 
!> client. After linking to libempire_api.{a|so} the Fortran client
!> communicates with EMPIRE via the module's public subroutines.
! 
!> \author Michael Andre
!> \author Altug Emiroglu
!> \date 07/11/2013
!> \comment Stefan Sicklinger The provided bind mechanism is part of 
!>                            FORTRAN 2003
!-----------------------------------------------------------------------
MODULE m_empire_api

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: EMPIRE_API_Connect, EMPIRE_API_sendMesh, &
            EMPIRE_API_sendDataField, EMPIRE_API_recvDataField, &
            EMPIRE_API_sendSignal_double, EMPIRE_API_recvSignal_double, &
            EMPIRE_API_recvConvergenceSignal, EMPIRE_API_printDataField, &
            EMPIRE_API_Disconnect, EMPIRE_API_getUserDefinedText!, writeDotMsh
            
  INTEGER, PARAMETER :: EMPIRE_API_NAME_STRING_LENGTH = 80
  
  
  INTERFACE

     ! -------------------------------------------------------------
     !> C-Binding: EMPIRE_API_Connect
     SUBROUTINE EMPIRE_API_Connect_c(emperor_input) &
       BIND(C, name="EMPIRE_API_Connect")
       USE iso_c_binding, ONLY : c_char
       CHARACTER(c_char), DIMENSION(*), INTENT(in) :: emperor_input !< input file
     END SUBROUTINE EMPIRE_API_Connect_c
     
     ! -------------------------------------------------------------
     !> C-Binding: EMPIRE_API_getUserDefinedText
     FUNCTION EMPIRE_API_getUserDefinedText_c(elementName) &
       BIND(C, name="EMPIRE_API_getUserDefinedText")
       USE iso_c_binding!, ONLY : c_char
       TYPE(c_ptr) EMPIRE_API_getUserDefinedText_c !< pointer to the XML element
       CHARACTER(c_char), DIMENSION(*), INTENT(in) :: elementName !< name of the XML element
     END FUNCTION EMPIRE_API_getUserDefinedText_c
       
     !--------------------------------------------------------------
     !> C-Binding: EMPIRE_API_sendMesh
     SUBROUTINE EMPIRE_API_sendMesh_c(mesh_name, num_nodes, num_elems, nodes, &
          node_ids, num_nodes_per_elem, elems) &
       BIND(C, name="EMPIRE_API_sendMesh")
       USE iso_c_binding, ONLY : c_char, c_int, c_double
       CHARACTER(c_char), DIMENSION(*), INTENT(in) :: mesh_name !< Mesh name.
       INTEGER(c_int), VALUE, INTENT(in) :: num_nodes !< Number of nodes.
       INTEGER(c_int), VALUE, INTENT(in) :: num_elems !< Number of elements.
       REAL(c_double), DIMENSION(*), INTENT(in) :: nodes !< Node coordinates.
       INTEGER(c_int), DIMENSION(*), INTENT(in) :: node_ids !< Node IDs.
       INTEGER(c_int), DIMENSION(*), INTENT(in) :: num_nodes_per_elem !< number of nodes at each element.
       INTEGER(c_int), DIMENSION(*), INTENT(in) :: elems !< Elements IDs.
     END SUBROUTINE EMPIRE_API_sendMesh_c

     !-------------------------------------------------------------
     !> C-Binding: EMPIRE_API_sendDataField
     SUBROUTINE EMPIRE_API_sendDataField_c(data_name, size, data_field) &
       BIND(C, name="EMPIRE_API_sendDataField")
       USE iso_c_binding, ONLY : c_char, c_int, c_double
       CHARACTER(c_char), DIMENSION(*), INTENT(in) :: data_name !< Data name.
       INTEGER(c_int), VALUE, INTENT(in) :: size !< Data array size.
       REAL(c_double), DIMENSION(*), INTENT(in) :: data_field !< Data array.
     END SUBROUTINE EMPIRE_API_sendDataField_c

     !-------------------------------------------------------------
     !> C-Binding: EMPIRE_API_recvDataField
     SUBROUTINE EMPIRE_API_recvDataField_c(data_name, size, data_field) &
       BIND(C, name="EMPIRE_API_recvDataField")
       USE iso_c_binding, ONLY : c_char, c_int, c_double
       CHARACTER(c_char), DIMENSION(*), INTENT(in) :: data_name !< Data name.
       INTEGER(c_int), VALUE, INTENT(in) :: size !< Data array size.
       REAL(c_double), DIMENSION(*), INTENT(out) :: data_field !< Data array.
     END SUBROUTINE EMPIRE_API_recvDataField_c

     !--------------------------------------------------------------
     !> C-Binding: EMPIRE_API_sendSignal_double
     SUBROUTINE EMPIRE_API_sendSignal_double_c(name, size, signal) &
       BIND(C, name="EMPIRE_API_sendSignal_double")
       USE iso_c_binding, ONLY : c_char, c_int, c_double
       CHARACTER(c_char), DIMENSION(*), INTENT(in) :: name !< Signal name.
       INTEGER(c_int), VALUE, INTENT(in) :: size !< Signal array size.
       REAL(c_double), DIMENSION(*), INTENT(in) :: signal !< Signal array.
     END SUBROUTINE EMPIRE_API_sendSignal_double_c

     !--------------------------------------------------------------
     !> C-Binding: EMPIRE_API_recvSignal_double
     SUBROUTINE EMPIRE_API_recvSignal_double_c(name, size, signal) &
       BIND(C, name="EMPIRE_API_recvSignal_double")
       USE iso_c_binding, ONLY : c_char, c_int, c_double
       CHARACTER(c_char), DIMENSION(*), INTENT(in) :: name !< Signal name.
       INTEGER(c_int), VALUE, INTENT(in) :: size !< Signal array size.
       REAL(c_double), DIMENSION(*), INTENT(out) :: signal !< Signal array.
     END SUBROUTINE EMPIRE_API_recvSignal_double_c

     !--------------------------------------------------------------
     !> C-Binding: EMPIRE_API_recvConvergenceSignal
     FUNCTION EMPIRE_API_recvConvergenceSignal_c() &
       BIND(C, name="EMPIRE_API_recvConvergenceSignal")
       USE iso_c_binding, ONLY : c_int
       INTEGER(c_int) EMPIRE_API_recvConvergenceSignal_c
     END FUNCTION EMPIRE_API_recvConvergenceSignal_c

     !--------------------------------------------------------------
     !> C-Binding: EMPIRE_API_printDataField
     SUBROUTINE EMPIRE_API_printDataField_c(name, size, data_field) &
       BIND(C, name="EMPIRE_API_printDataField")
       USE iso_c_binding, ONLY : c_char, c_int, c_double
       CHARACTER(c_char), DIMENSION(*), INTENT(in) :: name !< Data field name.
       INTEGER(c_int), VALUE, INTENT(in) :: size !< Data field size.
       REAL(c_double), DIMENSION(*), INTENT(in) :: data_field !< Data field array.
     END SUBROUTINE EMPIRE_API_printDataField_c

     !--------------------------------------------------------------
     !> C-Binding: EMPIRE_API_Disconnect
     SUBROUTINE EMPIRE_API_Disconnect_c() &
       BIND(C, name="EMPIRE_API_Disconnect")
     END SUBROUTINE EMPIRE_API_Disconnect_c

  END INTERFACE

CONTAINS

  ! ----------------------------------------------------------------
  !> Wrapper: EMPIRE_API_Connect
  SUBROUTINE EMPIRE_API_Connect(emperor_input)
    USE iso_c_binding, ONLY : c_char, c_null_char
    CHARACTER(len=*), INTENT(in) :: emperor_input
    CALL EMPIRE_API_Connect_c(emperor_input//c_null_char)
  END SUBROUTINE EMPIRE_API_Connect

  ! -------------------------------------------------------------
  !> Wrapper: EMPIRE_API_getUserDefinedText
  SUBROUTINE EMPIRE_API_getUserDefinedText(elementVal, elementName)
    USE iso_c_binding!, ONLY : c_char, c_null_char
    CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: elementVal !< value of the XML element
    CHARACTER(len=*), INTENT(in) :: elementName !< name of the XML element
    TYPE(c_ptr) :: c_p
    CHARACTER(c_char), POINTER :: f_p(:)
    INTEGER :: strLen

    c_p = EMPIRE_API_getUserDefinedText_c(elementName//c_null_char)
    call c_f_pointer(c_p,f_p,[EMPIRE_API_NAME_STRING_LENGTH])
    strLen = C_StrLen(f_p)
    allocate(character(len=strLen) :: elementVal)
    elementVal = transfer(f_p(1:strLen),elementVal)

  END SUBROUTINE EMPIRE_API_getUserDefinedText
  
  !-----------------------------------------------------------------
  !> Wrapper: EMPIRE_API_sendMesh
  SUBROUTINE EMPIRE_API_sendMesh(mesh_name, num_nodes, num_elems, nodes, &
       node_ids, num_nodes_per_elem, elems)
    USE iso_c_binding, ONLY : c_char, c_null_char
    CHARACTER(len=*), INTENT(in) :: mesh_name !< Mesh name.
    INTEGER, INTENT(in) :: num_nodes !< Number of nodes.
    INTEGER, INTENT(in) :: num_elems !< Number of elements.
    REAL(8), DIMENSION(:), INTENT(in) :: nodes !< Node coordinates.
    INTEGER, DIMENSION(:), INTENT(in) :: node_ids !< Node IDs.
    INTEGER, DIMENSION(:), INTENT(in) :: num_nodes_per_elem !< number of nodes at each element.
    INTEGER, DIMENSION(:), INTENT(in) :: elems !< Elements IDs.

    CALL EMPIRE_API_sendMesh_c(mesh_name//c_null_char, num_nodes, num_elems, nodes, &
         node_ids, num_nodes_per_elem, elems)
  END SUBROUTINE EMPIRE_API_sendMesh

  !-------------------------------------------------------------
  !> Wrapper: EMPIRE_API_sendDataField
  SUBROUTINE EMPIRE_API_sendDataField(data_name, size, data_field)
    USE iso_c_binding, ONLY : c_char, c_null_char
    CHARACTER(len=*), INTENT(in) :: data_name !< Data name.
    INTEGER, INTENT(in) :: size !< Data array size.
    REAL(8), DIMENSION(:), INTENT(in) :: data_field !< Data array.
    
    CALL EMPIRE_API_sendDataField_c(data_name//c_null_char, size, data_field)
  END SUBROUTINE EMPIRE_API_sendDataField
  
  !-------------------------------------------------------------
  !> Wrapper: EMPIRE_API_recvDataField
  SUBROUTINE EMPIRE_API_recvDataField(data_name, size, data_field)
    USE iso_c_binding, ONLY : c_char, c_null_char
    CHARACTER(len=*), INTENT(in) :: data_name !< Data name.
    INTEGER, VALUE, INTENT(in) :: size !< Data array size.
    REAL(8), DIMENSION(:), INTENT(out) :: data_field !< Data array.
    
    CALL EMPIRE_API_recvDataField_c(data_name//c_null_char, size, data_field)
  END SUBROUTINE EMPIRE_API_recvDataField
  
  !--------------------------------------------------------------
  !> Wrapper: EMPIRE_API_sendSignal_double
  SUBROUTINE EMPIRE_API_sendSignal_double(name, size, signal)
    USE iso_c_binding, ONLY : c_char, c_null_char
    CHARACTER(len=*), INTENT(in) :: name !< Signal name.
    INTEGER, INTENT(in) :: size !< Signal array size.
    REAL(8), DIMENSION(:), INTENT(in) :: signal !< Signal array.
    
    CALL EMPIRE_API_sendSignal_double_c(name//c_null_char, size, signal)
  END SUBROUTINE EMPIRE_API_sendSignal_double

  !--------------------------------------------------------------
  !> Wrapper: EMPIRE_API_recvSignal_double
  SUBROUTINE EMPIRE_API_recvSignal_double(name, size, signal)
    USE iso_c_binding, ONLY : c_char, c_null_char
    CHARACTER(len=*), INTENT(in) :: name !< Signal name.
    INTEGER, INTENT(in) :: size !< Signal array size.
    REAL(8), DIMENSION(:), INTENT(out) :: signal !< Signal array.
    
    CALL EMPIRE_API_recvSignal_double_c(name//c_null_char, size, signal)
  END SUBROUTINE EMPIRE_API_recvSignal_double

  !--------------------------------------------------------------
  !> Wrapper: EMPIRE_API_recvConvergenceSignal
  FUNCTION EMPIRE_API_recvConvergenceSignal()
    INTEGER EMPIRE_API_recvConvergenceSignal

    EMPIRE_API_recvConvergenceSignal = EMPIRE_API_recvConvergenceSignal_c()
    RETURN
  END FUNCTION EMPIRE_API_recvConvergenceSignal

  !--------------------------------------------------------------
  !> Wrapper: EMPIRE_API_printDataField
  SUBROUTINE EMPIRE_API_printDataField(name, size, data_field)
    USE iso_c_binding, ONLY : c_char, c_null_char
    CHARACTER(len=*), INTENT(in) :: name !< Data field name.
    INTEGER, INTENT(in) :: size !< Data field size.
    REAL(8), DIMENSION(:), INTENT(in) :: data_field !< Data field array.

    CALL EMPIRE_API_printDataField_c(name//c_null_char, size, data_field)
  END SUBROUTINE EMPIRE_API_printDataField
     
  !--------------------------------------------------------------
  !> Wrapper: EMPIRE_API_Disconnect
  SUBROUTINE EMPIRE_API_Disconnect()
    CALL EMPIRE_API_Disconnect_c()
  END SUBROUTINE EMPIRE_API_Disconnect
  
  !--------------------------------------------------------------
  !> Function: Returns the number of characters in a c_char array
  FUNCTION C_StrLen(cCharArray) RESULT(strLen)
    use iso_c_binding
    character(c_char), intent(in) :: cCharArray(:)
    integer :: strLen
    integer :: i

    DO i = 1, size(cCharArray)
      IF (cCharArray(i) == c_null_char) THEN
        strLen = i - 1
        RETURN
      END IF
    END DO
    strLen = i

  END FUNCTION C_StrLen

END MODULE m_empire_api
