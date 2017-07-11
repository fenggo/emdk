function torsion_angle(iatom,jatom,katom,latom)
   use mode, only: x
   implicit none

   integer      iatom, jatom,katom,latom
   real(8) ::   torsion_angle
   real(8)      r1(3),r2(3),r3(3),ra(3),rb(3)
   real(8)      radi,rsq1,rsq2
   real(8)      cos_ang
   real(8)   :: pi=3.1415926535897932

   ! compute torsion angle
   r1(1)=x(1,jatom)-x(1,iatom)
   r1(2)=x(2,jatom)-x(2,iatom)
   r1(3)=x(3,jatom)-x(3,iatom)

   r2(1)=x(1,katom)-x(1,jatom)
   r2(2)=x(2,katom)-x(2,jatom)
   r2(3)=x(3,katom)-x(3,jatom)

   r3(1)=x(1,latom)-x(1,katom)
   r3(2)=x(2,latom)-x(2,katom)
   r3(3)=x(3,latom)-x(3,katom)

   ! cross product of vectors
   ra(1) = r1(2)*r2(3) - r1(3)*r2(2)
   ra(2) = r1(3)*r2(1) - r1(1)*r2(3)
   ra(3) = r1(1)*r2(2) - r1(2)*r2(1)

   rb(1) = r2(2)*r3(3) - r2(3)*r3(2)
   rb(2) = r2(3)*r3(1) - r2(1)*r3(3)
   rb(3) = r2(1)*r3(2) - r2(2)*r3(1)

   !dot  product of vectors
   rsq1= ra(1)*ra(1)+ra(2)*ra(2)+ra(3)*ra(3)
   rsq2= rb(1)*rb(1)+rb(2)*rb(2)+rb(3)*rb(3)
   radi= sqrt(rsq1*rsq2)
   !print *,radi,iatom,jatom,katom,latom
   cos_ang=ra(1)*rb(1) + ra(2)*rb(2) + ra(3)*rb(3)
   !print *,cos_ang
   cos_ang=-cos_ang/radi
   !print *,cos_ang
   torsion_angle = 180.0 - acos(cos_ang)*180.0/pi
   !print *,dihe
   return
end function  torsion_angle
