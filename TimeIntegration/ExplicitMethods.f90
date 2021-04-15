module ExplicitMethods

    use Constants
    use PDEclass
    use ElementClass, only: Elem_t

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public subroutines
    public :: DGstepByRK3
    public :: DGstepByRK5

contains

!··············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 62.
!
!> @param[in]     tn   initial time of the integration step
!> @param[in]     dt   time step to integrate
!··············································································
subroutine DGstepByRK3(tn, dt)
    !* Arguments *!
    real(wp), intent(in) :: tn
    real(wp), intent(in) :: dt

    !* Local variables *!
    real(wp), save :: a(3) = [0.0_wp,       -5.0_wp/9.0_wp,  -153.0_wp/128.0_wp]
    real(wp), save :: b(3) = [0.0_wp,        1.0_wp/3.0_wp,   3.0_wp/4.0_wp    ]
    real(wp), save :: g(3) = [1.0_wp/3.0_wp, 15.0_wp/16.0_wp, 8.0_wp/15.0_wp   ]
    integer  :: i
    integer  :: j
    real(wp) :: t
    type(Elem_t), pointer :: elem

    ! Check that 'dt' /= 0
    if (dt <= 0.0_wp) return

    ! Step loop
    do i = 1, 3

        t = tn + b(i) * dt
        call PDE%globalTimeDerivative(t)

        ! Elements loop
        call PDE%mesh%elems%reset_last(0)
        do j = 1, PDE%mesh%elems%size()

            elem => PDE%mesh%elems%next()

            ! RK3 algorithm
            elem%G   = a(i) * elem%G + elem%PhiD
            elem%Phi = elem%Phi + g(i) * dt * elem%G

        end do

    end do

    nullify(elem)

end subroutine DGstepByRK3

!··············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Runge-Kutta, low storage method of fifth order.
!
!> @param[in]     tn   initial time of the integration step
!> @param[in]     dt   time step to integrate
!··············································································
subroutine DGstepByRK5(tn, dt)
    !* Arguments *!
    real(wp), intent(in) :: tn
    real(wp), intent(in) :: dt

    !* Local variables *!
    real(wp), save :: a(5) = [0.0_wp,             -0.4178904745_wp,   &
                              -1.192151694643_wp, -1.697784692471_wp, &
                              -1.514183444257_wp ]
    real(wp), save :: b(5) = [0.0_wp,             0.1496590219993_wp, &
                              0.3704009573644_wp, 0.6222557631345_wp, &
                              0.9582821306748_wp ]
    real(wp), save :: g(5) = [0.1496590219993_wp, 0.3792103129999_wp, &
                              0.8229550293869_wp, 0.6994504559488_wp, &
                              0.1530572479681_wp ]
    integer  :: i
    integer  :: j
    real(wp) :: t
    type(Elem_t), pointer :: elem

    ! Check that 'dt' /= 0
    if (dt <= 0.0_wp) return

    ! Step loop
    do i = 1, 5

        t = tn + b(i) * dt
        call PDE%globalTimeDerivative(t)

        ! Elements loop
        call PDE%mesh%elems%reset_last(0)
        do j = 1, PDE%mesh%elems%size()

            !---- Begin associate
            elem => PDE%mesh%elems%next()

            ! RK5 algorithm
            elem%G   = a(i) * elem%G + elem%PhiD
            elem%Phi = elem%Phi + g(i) * dt * elem%G

        end do

    end do

    nullify(elem)

end subroutine DGstepByRK5

end module ExplicitMethods
