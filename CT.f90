module modSystem

    implicit none
    character(len=300) fname,fnamechk
   
end module modSystem
!=======================================================================
program Main 

    call Charge_Transport
    
end program Main
!=======================================================================
subroutine Charge_Transport
    
    implicit none
    integer(kind=4) rec_char,Multi,Multi1,Multi2,i,j
    integer(kind=4) Nbas,Nbas1,Nbas2,Nocc,Nocc1,Nocc2
    real(kind=8),allocatable :: MoCu(:,:),En(:),Overlap(:,:),Hamilton(:,:)
    real(kind=8),allocatable :: MoCu0(:,:),MoCu1(:,:),En1(:),MoCu2(:,:),En2(:)
    character(len=10) t1
    character(len=8) d1
    character(len=50) sfname
    character(len=1) obt
   
    call Head_Prient(t1,d1)
    !-------------------------------------------------------------------
    write(*,*) "---------------------------------------------------"
    write(*,*) "Input for Dimer:"
    call Filename(sfname,rec_char,2)
    open(unit=100,file = sfname(1:rec_char)//"_CT.dat")
    write(100,*) "---------------------------------------------------"
    write(100,*) "Dimer:"
    call CountNbas(Nbas,Multi)
    allocate(MoCu(Nbas,Nbas),En(Nbas),Overlap(Nbas,Nbas), &
             Hamilton(Nbas,Nbas),MoCu0(Nbas,Nbas))
    obt = "O"
    call Find_Orbit(MoCu,En,Nbas,Nocc,obt)
    call Find_Overlap(Overlap,Nbas)
    call Calculate_Hamilton(Nbas,MoCu,En,Overlap,Hamilton)
    !-------------------------------------------------------------------
    write(*,*) "---------------------------------------------------"
    write(100,*) "---------------------------------------------------"
    write(*,*) "Input for Monomer1:"
    write(100,*) "Monomer1:"
    call Filename(sfname,rec_char,1)
    call CountNbas(Nbas1,Multi1)
    allocate(MoCu1(Nbas1,Nbas1),En1(Nbas1))
    obt = "O"
    call Find_Orbit(MoCu1,En1,Nbas1,Nocc1,obt)
    !-------------------------------------------------------------------
    write(*,*) "---------------------------------------------------"
    write(100,*) "---------------------------------------------------"
    write(*,*) "Input for Monomer2:"
    write(100,*) "Monomer2:"
    call Filename(sfname,rec_char,1)
    call CountNbas(Nbas2,Multi2)
    allocate(MoCu2(Nbas2,Nbas2),En2(Nbas2))
    obt = "O"
    call Find_Orbit(MoCu2,En2,Nbas2,Nocc2,obt)
    !-------------------------------------------------------------------
    if( Nbas1 + Nbas2 /= Nbas ) then
       write(*,*) "Stop : Error on Nbas, check the input!"
       stop
    end if
    !-------------------------------------------------------------------
    write(*,*) "---------------------------------------------------"
    write(100,*) "---------------------------------------------------"
    MoCu0(:,:) = 0.0d0
    do i = 1,Nbas1
       do j = 1,Nbas1
          MoCu0(j,i) = MoCu1(j,i)
       end do
    end do
    do i = Nbas1+1,Nbas
       do j = Nbas1+1,Nbas
          MoCu0(j,i) = MoCu2(j-Nbas1,i-Nbas1)
       end do
    end do
    call Calculate_Charge_Transport(Nbas,Nocc1,Nocc2,Nbas1,MoCu0,Overlap,Hamilton)
    !-------------------------------------------------------------------
    write(*,*)
    write(*,*) "Normal termination, the program have finished now !"
    write(100,*)
    write(100,*) "Normal termination, the program have finished now !"

end subroutine Charge_Transport
!=======================================================================
subroutine Calculate_Charge_Transport(Nbas,Nocc1,Nocc2,Nbas1,MoCu0, &
                                      Overlap,Fock)
    
    implicit none
    integer(kind=4) Nbas,Nocc1,Nocc2,Nbas1,i,j
    integer(kind=4) iStar,iEnd,jStar,jEnd,CTIout
    real(kind=8) MoCu0(Nbas,Nbas),Overlap(Nbas,Nbas)
    real(kind=8) Fock(Nbas,Nbas)
    real(kind=8) e1,e2,ee1,ee2,S12,J12,Je12
    real(kind=8) MTemp1(Nbas,Nbas),MTemp2(Nbas,Nbas)
    
    call Ao2Mo1(Nbas,MoCu0,Fock,MTemp1)
    call Ao2Mo1(Nbas,MoCu0,Overlap,MTemp2)
    
    write(*,*) " Press 1 to print important integrals."
    write(*,*) " Press 2 to print all integrals."
    read(*,*) CTIout
    
    write(*,*)
    write(*,*) "Calculating Charge Transfer Integral..."
    
    write(100,*)
    write(100,*) "================================================================"
    write(100,*)
    write(100,*) "       >>>>>>>>>>  Charge Transport Integral  <<<<<<<<<         "
    write(100,*)
    write(100,"(A,I5)") "   MO1 ----  HOMO-1 :",Nocc1-1
    write(100,"(A,I5)") "             HOMO   :",Nocc1                
    write(100,"(A,I5)") "             LUMO   :",Nocc1+1
    write(100,"(A,I5)") "             LUMO+1 :",Nocc1+2 
    write(100,*)
    write(100,"(A,I5)") "   MO2 ----  HOMO-1 :",Nocc2-1
    write(100,"(A,I5)") "             HOMO   :",Nocc2                
    write(100,"(A,I5)") "             LUMO   :",Nocc2+1
    write(100,"(A,I5)") "             LUMO+1 :",Nocc2+2
    write(100,*)
    write(100,*) "  MO1  MO2     ee1(eV)      ee2(eV)       Je(eV)    Je(kcal/mol)"
    write(100,*)
    
    if( CTIout == 1 ) then
       iStar = Nocc1-1
       iEnd  = Nocc1+2
       jStar = Nbas1+Nocc2-1
       jEnd  = Nbas1+Nocc2+2
    else
       iStar = 1
       iEnd  = Nbas1
       jStar = Nbas1+1
       jEnd  = Nbas
    end if
    
    do i = iStar,iEnd
       do j = jStar,jEnd
          e1 = MTemp1(i,i)
          e2 = MTemp1(j,j)
          J12 = MTemp1(j,i)
          S12 = MTemp2(j,i)
          ee1 = 0.5d0*((e1+e2)-2*J12*S12 + (e1-e2)*dsqrt(1-S12**2))/(1-S12**2)
          ee2 = 0.5d0*((e1+e2)-2*J12*S12 - (e1-e2)*dsqrt(1-S12**2))/(1-S12**2)
          Je12 = (J12 - 0.5d0*(e1 + e2)*S12)/(1 - S12**2)
          write(100,"(2I5,10F13.4)") i,j-Nbas1,ee1,ee2,Je12,Je12*23.05d0
       end do
    end do
		
end subroutine Calculate_Charge_Transport
!=======================================================================
subroutine Filename(sfname,rec_char,ijud)
    
    use modSystem,only : fname,fnamechk
    
    implicit none
    integer(kind=4) N,i,j,rec_char,ijud
    character(len=50) sfname
    logical file_exists
    
   if( ijud == 1 ) goto 30
        
20  write(*,*) " Please input the Gaussian log file :<e.g. out.log>"
    read(*,*) fname
    
    inquire(file = fname,exist = file_exists)   
    if(.not. file_exists) then
       write(*,*) ' Unable to open input file "',Trim(fname),'",input again.'
       write(*,*)
       goto 20
    end if
 
    N = len(fname)
    j = 0
    rec_char = 0
    
    do i = N,1,-1
       if( fname(i:i) == "/" ) then ! change "/" to "\" under windows.
          j = i
          exit
       end if
    end do

    do i = j+1,N
       if( fname(i:i) == "." ) exit
       rec_char = rec_char + 1
       sfname(rec_char:rec_char) = fname(i:i)
    end do
    
    fnamechk = sfname(1:rec_char)//".fchk"
    return

30  write(*,*) " Please input the Gaussian fchk file :<e.g. out.fchk>"
    read(*,*) fnamechk
    
    inquire(file = fnamechk,exist = file_exists)   
    if(.not. file_exists) then
       write(*,*) ' Unable to open input file "',Trim(fnamechk),'",input again.'
       write(*,*)
       goto 30
    end if
       
end subroutine Filename
!=======================================================================
subroutine CountNbas(Nbas,Multi)
    
    use modSystem,only : fnamechk
    
    implicit none
    integer(kind=4) Nbas,Multi,NeleA,NeleB
    character(len=20) ch1,ch(5)
  
    open(unit=10,file=fnamechk)
    
    do while( .true. ) 
       read(10,*) ch1
       if( ch1 == "Charge" ) exit
    end do
    read(10,*)
    read(10,*) 
    read(10,*) ch(1:5),NeleA
    read(10,*) ch(1:5),NeleB
    read(10,*) ch(1:5),Nbas
    
    Multi = NeleA-NeleB+1
    
    write(*,*)
    write(*,"(A,I5)") "  Spin Multiplicity:         ",Multi
    write(*,"(A,I5)") "  Number of basis functions: ",Nbas
    write(100,*)
    write(100,"(A,I5)") " Spin Multiplicity:         ",Multi
    write(100,"(A,I5)") " Number of basis functions: ",Nbas
    
    if( Multi > 1 ) then
       write(*,*) "Stop : Only Supports Close Shell !"
       stop
    end if
    
    close(10)
    
end subroutine CountNbas
!=======================================================================
subroutine Find_Orbit(MoCu,En,Nbas,Nocc,obt)
    
    use modSystem,only : fnamechk
    
    implicit none  
    integer(kind=4) Nbas,Nocc,NeleA,NeleB
    real(kind=8) MoCu(Nbas,Nbas),En(Nbas),Ent(Nbas)
    character(len=20) ch1,ch(5)
    character(len=1) obt
   
    open(unit=1101,file=fnamechk)
    write(*,*) " Readind Orbit..."
    !-------------------------------------------------------------------
    if( obt == "O" ) then
       do while( .true. ) 
          read(1101,*) ch1
          if( ch1 == "Alpha" ) exit
       end do
        
       read(1101,*) En
       read(1101,*)
       read(1101,*) MoCu    
    else if( obt == "A" ) then
       do while( .true. ) 
          read(1101,*) ch1
          if( ch1 == "Alpha" ) exit
       end do
        
       read(1101,*) En
       read(1101,*)
       read(1101,*) Ent
       read(1101,*)
       read(1101,*) MoCu
             
    else if( obt == "B" ) then
       read(1101,*) En
       read(1101,*)
       read(1101,*) MoCu
       read(1101,*)
       read(1101,*) MoCu    
    end if
    !-------------------------------------------------------------------
    rewind(1101)
    do while( .true. ) 
       read(1101,*) ch1
       if( ch1 == "Multiplicity " ) exit
    end do
    read(1101,*) 
    read(1101,*) ch(1:5),NeleA
    read(1101,*) ch(1:5),NeleB
    write(*,"(A,I5)") "  Number of alpha electrons: ",NeleA
    write(*,"(A,I5)") "  Number of beta  electrons: ",NeleB
    write(*,*)
    write(100,"(A,I5)") " Number of alpha electrons: ",NeleA
    write(100,"(A,I5)") " Number of beta  electrons: ",NeleB
    write(100,*)
    
    if( obt == "O" ) then
       Nocc = NeleA
    else if( obt == "A" ) then
       Nocc = NeleA
    else if( obt == "B" ) then
       Nocc = NeleB
    end if

    En(:) = 27.21138d0*En(:)
    
    close(1101)
    
end subroutine Find_Orbit
!=======================================================================
subroutine Find_Overlap(Overlap,Nbas)
    
    use modSystem,only : fname
    
    implicit none
    integer(kind=4) i,j,k,Nbas,iost
    real(kind=8) Overlap(Nbas,Nbas)
    character(len=15) ch1,ch2
   
    open(unit=12,file=fname)
    write(*,*) " Readind Overlap..."
    
    do while( .true. )
       read(12,*,iostat = iost) ch1,ch2  
       if(iost < 0) then 
           write(*,*)
           write(*,*) "Wrong input file, please check it ! "
           stop
       end if   
       if( ch2 == "Overlap" ) exit           
    end do

    do i = 0,Nbas/5-1
       read(12,*)
       read(12,*) k,Overlap(i*5+1,i*5+1)
       read(12,*) k,Overlap(i*5+2,i*5+1:i*5+2)
       read(12,*) k,Overlap(i*5+3,i*5+1:i*5+3)
       read(12,*) k,Overlap(i*5+4,i*5+1:i*5+4)
  
       do j = i*5 + 5,Nbas
          read(12,*) k,Overlap(j,i*5+1:i*5+5)
       end do
    end do
    
    read(12,*)
    do i = 1,Nbas-Nbas/5*5
       read(12,*) k,Overlap(Nbas/5*5+i,Nbas/5*5+1:Nbas/5*5+i)
    end do

    do i = 1,Nbas
       do j = i,Nbas
          Overlap(i,j) = Overlap(j,i)
       end do
    end do
    
    close(12)
    
end subroutine Find_Overlap
!=======================================================================
subroutine Calculate_Hamilton(Nbas,MoCu,En,Overlap,Hamilton)
     
    implicit none
    integer(kind=4) Nbas,n
    real(kind=8) Overlap(Nbas,Nbas),MoCu(Nbas,Nbas),En(Nbas)
    real(kind=8) EH(Nbas,Nbas),Hamilton(Nbas,Nbas)
    
    write(*,*) " Calculating Fock..."

    call Dgemm('T','N',Nbas,Nbas,Nbas,1.0d0,MoCu,Nbas,Overlap,Nbas,0.0d0,EH,Nbas)

    do n = 1,Nbas
       EH(n,1:Nbas) = EH(n,1:Nbas)*En(n)
    end do
    
    call MatrixSolve(MoCu,EH,Hamilton,Nbas,Nbas)
    
end subroutine Calculate_Hamilton
!=======================================================================
subroutine MatrixSolve(A,B,X,N,M)

   implicit none
   integer N,M,info,lda,ldb
   real(kind=8) A(N,N),B(N,M),X(N,M),ipiv(N),T(N,N)
   
   lda = N
   ldb = N
   
   T(:,:) = transpose(A(:,:))
   call dgesv(N,M,T,lda,ipiv,B,ldb,info)
   if(info /= 0) write(*,*) "Error in dgesv"
   X(:,:) = B(:,:)

end subroutine Matrixsolve
!=======================================================================
subroutine Ao2Mo1(Nbas,MoCu,Hmat,f_Mo)
     
    implicit none
    integer(kind=4) Nbas
    real(kind=8) MoCu(Nbas,Nbas),Hmat(Nbas,Nbas)
    real(kind=8) MTemp1(Nbas,Nbas),f_Mo(Nbas,Nbas)
    
    call dgemm('T','N',Nbas,Nbas,Nbas,1.0d0,MoCu,Nbas,Hmat,Nbas,0.0d0,MTemp1,Nbas)
    call dgemm('N','N',Nbas,Nbas,Nbas,1.0d0,MTemp1,Nbas,MoCu,Nbas,0.0d0,f_Mo,Nbas)
    
end subroutine Ao2Mo1
!=======================================================================
subroutine Head_Prient(t1,d1)

    implicit none
    character(len=10) t1
    character(len=8) d1
    
    call date_and_time( date = d1 , time = t1)
    
    write(*,*)
    write(*,*) "******************************************************* "
    write(*,*) "*   Charge Transport Integral Program for Gaussian    * "
    write(*,*) "*   This program need keywords 'nosymm iop(3/33=1)'   * "
    write(*,*) "*       Author : Warm_Cloud(3355196386@qq.com)        * "
    write(*,*) "******************************************************* "
    write(*,*)
    write(*,*) "                           ",d1(1:4),"-",d1(5:6),"-", &
                d1(7:8),"   ",t1(1:2),":",t1(3:4),":",t1(5:6)
    write(*,*)
    
end subroutine Head_Prient
!=======================================================================

