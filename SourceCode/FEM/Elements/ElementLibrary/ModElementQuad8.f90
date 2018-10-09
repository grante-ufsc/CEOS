!##################################################################################################
! This module has the attributes and methods of the four-node quadrilateral linear element with
! full integration
! ID -> Quad8 = 220
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications: Added Quad8 Element
! Date:         Author: Wagner Rupp
!##################################################################################################
module ElementQuad8

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use Element

	! Global variables within the module
	! -------------------------------------------------------------------------------------------
    real(8), pointer , dimension(:,:) :: NaturalCoordQuad8 => null()
    real(8), pointer , dimension(:)   :: WeightQuad8       => null()

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassElementQuad8: Attributes and methods of the element Quad8
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends(ClassElement) :: ClassElementQuad8

        ! Class Attributes: Inherited from ClassElement
        !--------------------------------------------------------------------------------------------
        contains
            ! Class Methods
            !--------------------------------------------------------------------------------------
            procedure :: GetProfile          => GetProfile_Quad8
            procedure :: GetGaussPoints      => GetGaussPoints_Quad8
            procedure :: GetNumberOfNodes    => GetNumberOfNodes_Quad8
            procedure :: GetShapeFunctions   => GetShapeFunctions_Quad8
            procedure :: GetDifShapeFunctions=> GetDifShapeFunctions_Quad8
            procedure :: AllocateGaussPoints => AllocateGaussPointsParameters_Quad8
            procedure :: IntegrateLine       => IntegrateLine_Quad8

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        subroutine GetProfile_Quad8(this,Profile)
            class(ClassElementQuad8)::this
            type(ClassElementProfile)::Profile

            call Profile % LoadProfile( 0                 , &
            NumberOfNodes = 8                              , &
            IsQuadratic = .true.                          , &
            GeometryType = GeometryTypes % Quadrilateral   , &
            FullIntegrationCapable = .true.                , &
            MeanDilatationCapable= .false.                  , &
            ElementDimension = 2 )

        end subroutine

        !==========================================================================================
        ! Method GetGaussPoints_Quad8:  This method points to the natural coordinates and weights
        ! used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetGaussPoints_Quad8(this, NaturalCoord, Weight)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementQuad8) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer, dimension(:,:)  :: NaturalCoord
            real(8) , pointer, dimension(:)    :: Weight

		    !************************************************************************************

		    !************************************************************************************
            ! POINT TO QUAD8 METHODS
		    !************************************************************************************

            NaturalCoord => NaturalCoordQuad8
            Weight       => WeightQuad8

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetNumberOfNodes_Quad8:  This method returns the number of nodes of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        function GetNumberOfNodes_Quad8(this) result(nNodes)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementQuad8) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer  :: nNodes

		    !************************************************************************************

		    !************************************************************************************
            ! NUMBER OF NODES - QUAD8
		    !************************************************************************************

            nNodes = 8

		    !************************************************************************************

        end function
        !==========================================================================================

        !==========================================================================================
        ! Method GetShapeFunctions_Quad8:  This method returns the shape funtions of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetShapeFunctions_Quad8(this , NaturalCoord , ShapeFunctions )

			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementQuad8) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(inout) :: ShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8)             :: xi , eta
            real(8) , parameter :: R1=1.0d0

	!************************************************************************************--

      	!************************************************************************************
            ! SHAPE FUNTIONS - QUAD8
      	!************************************************************************************

            xi = NaturalCoord(1)
            eta = NaturalCoord(2)

            ShapeFunctions(1) = - 1.0d0 / 4.0d0 * (R1 - xi) * (R1 + eta) * (+R1 + xi - eta)
            ShapeFunctions(2) = - 1.0d0 / 4.0d0 * (R1 + xi) * (R1 + eta) * (+R1 - xi - eta)            
            ShapeFunctions(3) = + 1.0d0 / 4.0d0 * (R1 + xi) * (R1 - eta) * (-R1 + xi - eta)
            ShapeFunctions(4) = + 1.0d0 / 4.0d0 * (R1 - xi) * (R1 - eta) * (-R1 - xi - eta)
            
            ShapeFunctions(5) = + 1.0d0 / 2.0d0 * ( R1 + xi ) * ( R1 - xi  ) * ( R1 + eta )                 
            ShapeFunctions(6) = + 1.0d0 / 2.0d0 * ( R1 + xi ) * ( R1 + eta ) * ( R1 - eta )            
            ShapeFunctions(7) = + 1.0d0 / 2.0d0 * ( R1 - xi ) * ( R1 + xi  ) * ( R1 - eta )
            ShapeFunctions(8) = + 1.0d0 / 2.0d0 * ( R1 - xi ) * ( R1 - eta ) * ( R1 + eta )

      	    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetDifShapeFunctions_Quad8:  This method returns the shape funtions derivatives
        ! of the element.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetDifShapeFunctions_Quad8(this , NaturalCoord , DifShapeFunctions )

			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementQuad8) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: DifShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8)             :: xi , eta
            real(8) , parameter :: R1 = 1.0d0 , R2 = 1.0d0/2.0d0 , R4 = 1.0d0/4.0d0

            !************************************************************************************
            xi = NaturalCoord(1) ; eta = NaturalCoord(2)

      	!************************************************************************************
            ! SHAPE FUNTIONS DERIVATIVE- QUAD8
      	!************************************************************************************

            DifShapeFunctions(1,1) =-((eta - 2.0d0*xi)*(eta + 1.0d0))/4.0d0 ;
            DifShapeFunctions(1,2) =-((xi - 1.0d0)*(2.0d0*eta - xi))/4.0d0 ;

            DifShapeFunctions(2,1) =((eta + 2.0d0*xi)*(eta + 1.0d0))/4.0d0;
            DifShapeFunctions(2,2) = ((2.0d0*eta + xi)*(xi + 1.0d0))/4.0d0;
 
            DifShapeFunctions(3,1) =((eta - 2.0d0*xi)*(eta - 1.0d0))/4.0d0; 
            DifShapeFunctions(3,2) = ((xi + 1.0d0)*(2.0d0*eta - xi))/4.0d0;
 
            DifShapeFunctions(4,1) = -((eta + 2.0d0*xi)*(eta - 1.0d0))/4.0d0; 
            DifShapeFunctions(4,2) = -((2.0d0*eta + xi)*(xi - 1.0d0))/4.0d0;
 
            DifShapeFunctions(5,1) =  -xi*(eta + 1.0d0); 
            DifShapeFunctions(5,2) =  1.0d0/2.0d0 - xi ** 2.0d0/2.0d0;
 
            DifShapeFunctions(6,1) = -((eta - 1.0d0)*(eta + 1.0d0))/2.0d0; 
            DifShapeFunctions(6,2) = -eta*(xi + 1.0d0);
 
            DifShapeFunctions(7,1) = xi*(eta - 1.0d0); 
            DifShapeFunctions(7,2) = xi ** 2.0d0 /2.0d0 - 1.0d0/2.0d0;
 
            DifShapeFunctions(8,1) = ((eta - 1.0d0)*(eta + 1.0d0))/2.0d0;
            DifShapeFunctions(8,2) = eta*(xi - 1.0d0);

        !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method AllocateGaussPointsParameters_Quad8: This method returns the natural coordinates
        ! and weights used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine AllocateGaussPointsParameters_Quad8(this,nGP)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementQuad8) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(inout) :: nGP

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: x

		    !************************************************************************************

		    !************************************************************************************
            ! PARAMETERS OF GAUSS POINTS - QUAD8
		    !************************************************************************************

            !Number of Gauss Points
            nGP=9

            if (associated(NaturalCoordQuad8)) return
            allocate( NaturalCoordQuad8(nGP,2) , WeightQuad8(nGP) )

            x = dsqrt(0.6d0)

            NaturalCoordQuad8(1,:)=[-x,+x]
            NaturalCoordQuad8(2,:)=[0.0d0, +x]
            NaturalCoordQuad8(3,:)=[+x, +x]

            NaturalCoordQuad8(4,:)=[-x, 0.0d0]
            NaturalCoordQuad8(5,:)=[0.0d0, 0.0d0]
            NaturalCoordQuad8(6,:)=[+x, 0.0d0]

            NaturalCoordQuad8(7,:)=[-x, -x]
            NaturalCoordQuad8(8,:)=[0.0d0, -x]
            NaturalCoordQuad8(9,:)=[+x, -x]

            WeightQuad8(1)= 25.0d0 / 81.0d0 
            WeightQuad8(3)= 25.0d0 / 81.0d0 
            WeightQuad8(7)= 25.0d0 / 81.0d0 
            WeightQuad8(9)= 25.0d0 / 81.0d0 
            
            WeightQuad8(2)= 40.0d0 / 81.0d0 
            WeightQuad8(4)= 40.0d0 / 81.0d0 
            WeightQuad8(6)= 40.0d0 / 81.0d0 
            WeightQuad8(8)= 40.0d0 / 81.0d0 
            
            WeightQuad8(5)= 64.0d0 / 81.0d0 
            

		    !************************************************************************************

            end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine IntegrateLine_Quad8(this,LineNodes,t,F)
            use MathRoutines
            use Nodes
            implicit none
            class(ClassElementQuad8):: this
            type(ClassElementNodes) , dimension(:) :: LineNodes
            real(8) , dimension(:) :: F
            real(8)  :: t

            real(8) :: L
            real(8) :: Ftotal

            print *, "error"
        end subroutine
        !==========================================================================================

end module
