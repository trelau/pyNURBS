module icurve
use math, only: cross
use evaluate, only: curve_point, rat_surface_derivs
use intersect_surface, only: refine_spi_point, refine_ssi_point
implicit none

private
public :: icurve_eval

contains

subroutine icurve_eval(u, n1, p1, uk1, cpw1, n2, p2, uk2, cpw2, &
                       sn1, sp1, suk1, sm1, sq1, svk1, spw1, &
                       sn2, sp2, suk2, sm2, sq2, svk2, spw2, &
                       p0, pnorm, vx, vy, is_planar, tol, pnt, cu)
    !> Evaluate a point on an intersection curve.

    !f2py intent(in) u, n1, p1 uk1, cpw1, n2, p2, uk2, cpw2
    !f2py intent(in) sn1, sp1, suk1, sm1, sq1, svk1, spw1
    !f2py intent(in) sn2, sp2, suk2, sm2, sq2, svk2, spw2
    !f2py intent(in) p0, pnorm
    !f2py intent(in) is_planar, tol
    
    !f2py intent(out) pnt, cu
    
    !f2py depend(n1) cpw1
    !f2py depend(n1, p1) uk1
    !f2py depend(n2) cpw2
    !f2py depend(n2, p2) uk2
    
    !f2py depend(sn1, sm1) spw1
    !f2py depend(sn1, sp1) suk1
    !f2py depend(sm1, sq1) svk1
    
    !f2py depend(sn2, sm2) spw2
    !f2py depend(sn2, sp2) suk2
    !f2py depend(sm2, sq2) svk2
    
    ! Required input
    logical, intent(in) :: is_planar
    integer, intent(in) :: n1, p1, n2, p2, sn1, sp1, sm1, sq1
    integer, intent(in) :: sn2, sp2, sm2, sq2
    double precision, intent(in) :: u, tol
    double precision, intent(in) :: uk1(0:n1 + p1 + 1)
    double precision, intent(in) :: cpw1(0:n1, 0:3)
    double precision, intent(in) :: uk2(0:n2 + p2 + 1)
    double precision, intent(in) :: cpw2(0:n2, 0:3)
    
    double precision, intent(in) :: suk1(0:sn1 + sp1 + 1)
    double precision, intent(in) :: svk1(0:sm1 + sq1 + 1)
    double precision, intent(in) :: spw1(0:sn1, 0:sm1, 0:3)
    
    double precision, intent(in) :: suk2(0:sn2 + sp1 + 1)
    double precision, intent(in) :: svk2(0:sm2 + sq2 + 1)
    double precision, intent(in) :: spw2(0:sn2, 0:sm2, 0:3)
    double precision, intent(in) :: p0(0:2), pnorm(0:2), vx(0:2), vy(0:2)
    
    ! Output
    double precision, intent(out) :: pnt(0:2), cu(0:2)
    
    ! Working
    double precision :: u1, v1, u2, v2
    double precision :: uv1(0:2), uv2(0:2)
    double precision :: pi(0:2)
    double precision :: skl(0:2, 0:2, 0:2)
    double precision :: du1(0:2), dv1(0:2), du2(0:2), dv2(0:2)
    double precision :: vn1(0:2), vn2(0:2)
    



    ! Evaluate 2-D curves.
    call curve_point(n1, p1, uk1, cpw1, u, uv1)
    call curve_point(n2, p2, uk2, cpw2, u, uv2)
    
    ! Refine points.
    if (is_planar) then
        call refine_spi_point(sn1, sp1, suk1, sm1, sq1, svk1, spw1, &
                              p0, pnorm, vx, vy, uv1(0), uv1(1), tol, &
                              u1, v1, u2, v2, pi)
    else
        call refine_ssi_point(n1, sp1, suk1, sm1, sq1, svk1, spw1, &
                              n2, sp2, suk2, sm2, sq2, svk2, spw2, &
                              uv1(0), uv1(1), uv2(0), uv2(1), tol, &
                              u1, v1, u2, v2, pi)
    end if
    
    ! Evaluate point on surface 1 and first derivatives.
    call rat_surface_derivs(sn1, sp1, suk1, sm1, sq1, svk1, spw1, &
                            u1, v1, 2, skl)
    pnt = skl(0, 0, :)
    du1 = skl(1, 0, :)
    dv1 = skl(0, 1, :)
    ! Cross product to get surface norm.
    vn1 = cross(du1, dv1)
    if (is_planar) then
        vn2 = pnorm
    else
        call rat_surface_derivs(sn2, sp2, suk2, sm2, sq2, svk2, spw2, &
                                u2, v2, 2, skl)
        du2 = skl(1, 0, :)
        dv2 = skl(0, 1, :)
        vn2 = cross(du2, dv2)
    end if

    ! Cross product of normal vectors to get derivative.
    cu = cross(vn1, vn2)

end subroutine icurve_eval

end module icurve