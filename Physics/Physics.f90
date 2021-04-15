!*******************************************************************************
!  MODULE: Physics
!
!> @author
!> Andres Mateo
!
!> @brief
!> Wraps all the physics-related modules.
!*******************************************************************************

module Physics

    ! Only the compressible Navier-Stokes equations are implemented
    use PhysicsStorage_NS
    use Physics_NS
    use RiemannSolvers_NS
    use ArtificialViscosity_NS
    use ManufacturedSolutions_NS

end module Physics
