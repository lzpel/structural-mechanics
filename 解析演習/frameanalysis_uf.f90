subroutine elemstiff(nelem,nnode,pos,ine,ielem,EA,EI,se)
  implicit none
  integer nelem,nnode
  real(8) pos(2,nnode),EA,EI,se(6,6)
  integer ine(2,nelem),inode,jnode,ielem,i,j,k
  real(8) lngth,dx(2),cs,sn,r(6,6),sr(6,6)
  inode=ine(1,ielem)
  jnode=ine(2,ielem)
  dx(1)=pos(1,jnode)-pos(1,inode)
  dx(2)=pos(2,jnode)-pos(2,inode)
  lngth=sqrt(dx(1)**2+dx(2)**2)
!!!***Compute[SM]***
!!$  se(1,1)=???
!!$  ???
  
  do i=1,6
     do j=1,i-1
        se(i,j)=se(j,i)
     end do
  end do
!!!***Compute[R]***
  cs=dx(1)/lngth
  sn=dx(2)/lngth
  r(:,:)=0.d0
!!$  r(1,1)=???
!!$  ???
!!!S*R
!!$  sr(:,:)=0.d0
!!$  do i=1,6
!!$     do j=1,6
!!$        do k=1,6
!!$           sr(i,j)=???
!!$        end do
!!$     end do
!!$  end do
!!!R^t*S*R
!!$  ???
end subroutine elemstiff

subroutine globstiff(ielem,ine,s,se,nelem,nnode,mapping)
  integer nelem,ielem, ine(2,nelem),nnode,mapping(3,nnode)
  real(8) se(6,6),s(3*nnode,3*nnode)
  integer:: kk(6)
  integer:: i, j, k, l
  
  inode=ine(1,ielem)
  jnode=ine(2,ielem)
!!$  kk(1)=???
!!$  ???
  
!!$  do i=1, 6
!!$     k=???
!!$     do j=1, 6
!!$        l=???
!!$        s(k,l)=s(k,l)+se(i,j)
!!$     end do
!!$  end do
end subroutine globstiff

subroutine setbc(nnode,ibc,s,force,u,mapping)
  implicit none
  integer nnode,ibc(3,nnode),mapping(3,nnode)
  real(8) force(3,nnode),s(3*nnode,3*nnode),u(3*nnode)
  integer:: i, j,inode,k

!!$  do inode=1,nnode
!!$     do i=1,3
!!$        k=???
!!$        u(k)=???
!!$     end do
!!$  end do
!!$  
!!$  do inode=1,nnode
!!$     do i=1,3
!!$        if(ibc(i,inode)==1) then
!!$           k=???
!!$           do j=1,3*nnode
!!$              s(k,j)=???
!!$              s(j,k)=???
!!$           end do
!!$           s(k,k)=???
!!$           u(k)=???
!!$        end if
!!$     end do
!!$  end do

end subroutine setbc


subroutine solvele(a,x,n)
  implicit none
  integer n
  real(8) a(n,n),x(n)
  integer i,j,k
  
  do i=1,n-1
     do j=i+1,n
        x(j)=x(j)-x(i)*a(j,i)/a(i,i)
        do k=i+1,n
           a(j,k)=a(j,k)-a(i,k)*a(j,i)/a(i,i)
        end do
     end do
  end do

  x(n)=x(n)/a(n,n)
  do i=n-1,1,-1
     do j=i+1,n
        x(i)=x(i)-a(i,j)*x(j)
     end do
     x(i)=x(i)/a(i,i)
  end do
end subroutine solvele

subroutine compstrain(nelem,ine,pos,disp,hght,EI,nnode)
  implicit none
  integer nelem,nnode,ine(2,nelem)
  real(8) pos(2,nnode),disp(3,nnode),hght(nelem),EI(nelem)
  real(8) lngth,dx(2),cs,sn
  real(8) v(2),psi(2),c2,c3,bm,tops
  integer inode,jnode,ielem
  do ielem=1,nelem
     inode=ine(1,ielem)
     jnode=ine(2,ielem)
     dx(1)=pos(1,jnode)-pos(1,inode)       
     dx(2)=pos(2,jnode)-pos(2,inode)
     lngth=sqrt(dx(1)**2+dx(2)**2)
     cs=dx(1)/lngth
     sn=dx(2)/lngth
     v(1)=-sn*disp(1,inode)+cs*disp(2,inode)
     v(2)=-sn*disp(1,jnode)+cs*disp(2,jnode)
     psi(1)=disp(3,inode)
     psi(2)=disp(3,jnode)
!!$     c2=???
!!$     c3=???
!!$     bm=???
!!$     tops=???
!!$     write(*,*) ielem,bm,tops
  end do
end subroutine compstrain

subroutine outputd(nelem,nnode,ine,pos,disp)
  implicit none
  integer nelem,nnode, ine(2,nelem)
  real(8) pos(2,nnode), disp(3,nnode)
  real(8) lngth, dx(2), cs, sn
  real(8) a0, a1, c0, c1, c2, c3
  real(8) u(2), v(2), psi(2)
  real(8) xi, xini, yini, ubar, vbar, udef, vdef
  integer j, inode, jnode, ielem
  open(200,file='deformation.dat') !!@\label{openwrite}@

  do ielem = 1, nelem
     inode=ine(1,ielem)
     jnode=ine(2,ielem)
     dx(1)=pos(1,jnode)-pos(1,inode)
     dx(2)=pos(2,jnode)-pos(2,inode)
     lngth=sqrt(dx(1)**2+dx(2)**2)
     cs=dx(1)/lngth
     sn=dx(2)/lngth

     u(1)=cs*disp(1,inode)+sn*disp(2,inode)
     v(1)=-sn*disp(1,inode)+cs*disp(2,inode)
     psi(1)=disp(3,inode)
     u(2)=cs*disp(1,jnode)+sn*disp(2,jnode)
     v(2)=-sn*disp(1,jnode)+cs*disp(2,jnode)
     psi(2)=disp(3,jnode)

     a0=u(1)
     a1=-u(1)/lngth+u(2)/lngth
     c0=v(1)
     c1=psi(1)
     c2=-3.d0*v(1)/lngth**2-2.d0*psi(1)/lngth &
          & +3.d0*v(2)/lngth**2-psi(2)/lngth
     c3=2.d0*v(1)/lngth**3+psi(1)/lngth**2 &
          & -2.d0*v(2)/lngth**3+psi(2)/lngth**2
     do j=1,11
        xi=lngth/10*(j-1)
        xini=dx(1)/10*(j-1)+pos(1,inode)
        yini=dx(2)/10*(j-1)+pos(2,inode)
        ubar=a0+a1*xi
        vbar=c0+(c1+(c2+c3*xi)*xi)*xi
!!$        udef=???
!!$        vdef=???
        write(200,*) xini, yini, udef, vdef
     end do
     write(200,*)
  end do
  close(200)
end subroutine outputd

program frame_analysis
  implicit none
  integer nnode,nelem
  integer,allocatable::ine(:,:),ibc(:,:),mapping(:,:)
  real(8),allocatable::pos(:,:),force(:,:),disp(:,:)
  real(8),allocatable::hght(:),wdth(:),ym(:),EA(:),EI(:)
  real(8),allocatable::s(:,:),u(:)
  real(8) se(6,6)
  integer ielem,inode,i,j
!!!===INPUTDATA===
  open(100,file='kunoji.dat')
  read(100,*)nnode,nelem
  allocate(pos(2,nnode),force(3,nnode),disp(3,nnode),mapping(3,nnode))
  allocate(hght(nelem),wdth(nelem),ym(nelem),EA(nelem),EI(nelem))
  allocate(s(3*nnode,3*nnode),u(3*nnode))
  allocate(ine(2,nelem),ibc(3,nnode))

  do inode=1,nnode
     read(100,*)(pos(i,inode),i=1,2)
  end do
  
  j=0
  do inode=1,nnode
     do i=1,3
!!$     ???
        mapping(i,inode)=j
     end do
  end do

  do ielem=1,nelem
     read(100,*)(ine(i,ielem),i=1,2)
  end do
  
!!$  do ielem=1,nelem
!!$     read(100,*)???
!!$     EA(ielem)=???
!!$     EI(ielem)=???
!!$  end do
!!$
!!$  do inode=1,nnode
!!$     read(100,*)???
!!$  end do
!!$
!!$  do inode=1,nnode
!!$     read(100,*)???
!!$  end do
!!$  close(100)
!!$  
!!$  !!!=== COMPUTE STIFFNESS MATRIX ===
!!$  s(:,:)=0.d0
!!$  do ielem=1, nelem
!!$     call elemstiff(nelem,nnode,pos,ine,ielem,EA(ielem),EI(ielem),se)
!!$     call globstiff(ielem,ine,s,se,nelem,nnode,mapping)
!!$  end do
!!$  do i=1,9
!!$     write(*,'(9e10.3e1)') (s(i,j),j=1,9)
!!$  end do
!!$  write(*,*)
!!$  
!!$!!!=== SET BOUNDARY CONDITIONS ===
!!$  call setbc(nnode,ibc,s,force,u,mapping)
!!$  do i=1,9
!!$     write(*,'(9e10.3e1)') (s(i,j),j=1,9)
!!$  end do
!!$
!!$!!!=== SOLVE LINEAR EQUATIONS ===
!!$  ???
!!$  do inode=1,nnode
!!$     do i=1,3
!!$        disp(i,inode)=???
!!$     end do
!!$  end do
!!$
!!$!!!=== COMPUTE STRAINS AT TOP(OUTSIDE) OF THE BEAMS ===
!!$  call compstrain(nelem,ine,pos,disp,hght,ei,nnode)
!!$
!!$!!!=== OUTPUT DEFORMATION ===
!!$  call outputd(nelem,nnode,ine,pos,disp)

end program frame_analysis
