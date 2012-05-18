      subroutine ehg182(i)
      integer i
      if(i.eq.100) print *,'wrong version number in lowesd.  Probably ty
     +po in caller.'
      if(i.eq.101) print *,'d>dMAX in ehg131.  Need to recompile with in
     +creased dimensions.'
      if(i.eq.102) print *,'liv too small.   (Discovered by lowesd)'
      if(i.eq.103) print *,'lv too small.    (Discovered by lowesd)'
      if(i.eq.104) print *,'alpha too small.  fewer data values than deg
     +rees of freedom.'
      if(i.eq.105) print *,'k>d2MAX in ehg136.  Need to recompile with i
     +ncreased dimensions.'
      if(i.eq.106) print *,'lwork too small'
      if(i.eq.107) print *,'invalid value for kernel'
      if(i.eq.108) print *,'invalid value for ideg'
      if(i.eq.109) print *,'lowstt only applies when kernel=1.'
      if(i.eq.110) print *,'not enough extra workspace for robustness ca
     +lculation'
      if(i.eq.120) print *,'zero-width neighborhood. make alpha bigger'
      if(i.eq.121) print *,'all data on boundary of neighborhood. make a
     +lpha bigger'
      if(i.eq.122) print *,'extrapolation not allowed with blending'
      if(i.eq.123) print *,'ihat=1 (diag L) in l2fit only makes sense if
     + z=x (eval=data).'
      if(i.eq.171) print *,'lowesd must be called first.'
      if(i.eq.172) print *,'lowesf must not come between lowesb and lowe
     +se, lowesr, or lowesl.'
      if(i.eq.173) print *,'lowesb must come before lowese, lowesr, or l
     +owesl.'
      if(i.eq.174) print *,'lowesb need not be called twice.'
      if(i.eq.180) print *,'nv>nvmax in cpvert.'
      if(i.eq.181) print *,'nt>20 in eval.'
      if(i.eq.182) print *,'svddc failed in l2fit.'
      if(i.eq.183) print *,'didnt find edge in vleaf.'
      if(i.eq.184) print *,'zero-width cell found in vleaf.'
      if(i.eq.185) print *,'trouble descending to leaf in vleaf.'
      if(i.eq.186) print *,'insufficient workspace for lowesf.'
      if(i.eq.187) print *,'insufficient stack space'
      if(i.eq.188) print *,'lv too small for computing explicit L'
      if(i.eq.191) print *,'computed trace L was negative; something is 
     +wrong!'
      if(i.eq.192) print *,'computed delta was negative; something is wr
     +ong!'
      if(i.eq.193) print *,'workspace in loread appears to be corrupted'
      if(i.eq.194) print *,'trouble in l2fit/l2tr'
      if(i.eq.195) print *,'only constant, linear, or quadratic local mo
     +dels allowed'
      if(i.eq.196) print *,'degree must be at least 1 for vertex influen
     +ce matrix'
      if(i.eq.999) print *,'not yet implemented'
      print *,'Assert failed, error code ',i
      stop
      end
      subroutine ehg183(s,i,n,inc)
      character*(*) s
      integer n, inc, i(inc,n), j
      print *,s,(i(1,j),j=1,n)
      end
      subroutine ehg184(s,x,n,inc)
      character*(*) s
      integer n, inc, j
      double precision x(inc,n)
      print *,s,(x(1,j),j=1,n)
      end
      subroutine losave(iunit,iv,liv,lv,v)
      integer execnt,iunit,liv,lv
      integer iv(liv)
      DOUBLE PRECISION v(lv)
      external ehg167
      save execnt
      data execnt /0/
      execnt=execnt+1
      call ehg167(iunit,iv(2),iv(4),iv(5),iv(6),iv(14),v(iv(11)),iv(iv(7
     +)),v(iv(12)),v(iv(13)))
      return
      end
      subroutine ehg167(iunit,d,vc,nc,nv,nvmax,v,a,xi,vval)
      integer iunit,d,vc,nc,nv,a(nc),magic,i,j
      DOUBLE PRECISION v(nvmax,d),xi(nc),vval(0:d,nv)
      write(iunit,*)d,nc,nv
      do 10 i=1,d
10      write(iunit,*)v(1,i),v(vc,i)
      j = 0
      do 20 i=1,nc
        if(a(i).ne.0)then
          write(iunit,*)a(i),xi(i)
        else
          write(iunit,*)a(i),j
        end if
20    continue
      do 30 i=1,nv
30      write(iunit,*)(vval(j,i),j=0,d)
      end
      subroutine lohead(iunit,d,vc,nc,nv)
      integer iunit,d,vc,nc,nv
      read(iunit,*)d,nc,nv
      vc = 2**d
      end
      subroutine loread(iunit,d,vc,nc,nv,iv,liv,lv,v)
      integer bound,d,execnt,iunit,liv,lv,nc,nv,vc
      integer iv(liv)
      DOUBLE PRECISION v(lv)
      external ehg168,ehg169,ehg182
      save execnt
      data execnt /0/
      execnt=execnt+1
      iv(28)=173
      iv(2)=d
      iv(4)=vc
      iv(14)=nv
      iv(17)=nc
      iv(7)=50
      iv(8)=iv(7)+nc
      iv(9)=iv(8)+vc*nc
      iv(10)=iv(9)+nc
      bound=iv(10)+nc
      if(.not.(bound-1.le.liv))then
         call ehg182(102)
      end if
      iv(11)=50
      iv(13)=iv(11)+nv*d
      iv(12)=iv(13)+(d+1)*nv
      bound=iv(12)+nc
      if(.not.(bound-1.le.lv))then
         call ehg182(103)
      end if
      call ehg168(iunit,d,vc,nc,nv,nv,v(iv(11)),iv(iv(7)),v(iv(12)),v(iv
     +(13)))
      call ehg169(d,vc,nc,nc,nv,nv,v(iv(11)),iv(iv(7)),v(iv(12)),iv(iv(8
     +)),iv(iv(9)),iv(iv(10)))
      return
      end
      subroutine ehg168(iunit,d,vc,nc,nv,nvmax,v,a,xi,vval)
      integer iunit,d,vc,nc,nv,a(nc),magic,i,j
      DOUBLE PRECISION v(nvmax,d),xi(nc),vval(0:d,nv)
      do 10 i=1,d
10      read(iunit,*)v(1,i),v(vc,i)
      do 20 i=1,nc
20      read(iunit,*)a(i),xi(i)
      do 30 i=1,nv
30      read(iunit,*)(vval(j,i),j=0,d)
      end
      subroutine ehg170(k,d,vc,nv,nvmax,nc,ncmax,a,c,hi,lo,v,vval,xi)
      integer d,execnt,i,j,nc,ncmax,nv,nvmax,vc
      integer a(ncmax),c(vc,ncmax),hi(ncmax),lo(ncmax)
      double precision v(nvmax,d),vval(0:d,nvmax),xi(ncmax)
      save execnt
      data execnt /0/
      execnt=execnt+1
      write(k,*)'      double precision function loeval(z)'
      write(k,50)d
      write(k,*)'      integer d,vc,nv,nc'
      write(k,51)nc,vc,nc
      write(k,52)nc,nc
      write(k,53)nv,d
      write(k,54)d,nv
      write(k,55)nc
      write(k,56)
      write(k,57)d,vc,nv,nc
50    format('      double precision z(',i2,')')
51    format('      integer a(',i5,'), c(',i3,',',i5,')')
52    format('      integer hi(',i5,'), lo(',i5,')')
53    format('      double precision v(',i5,',',i2,')')
54    format('      double precision vval(0:',i2,',',i5,')')
55    format('      double precision xi(',i5,')')
56    format('      double precision ehg128')
57    format('      data d,vc,nv,nc /',i2,',',i3,',',i5,',',i5,'/')
      do 3 i=1,nc
         write(k,58)i,a(i)
58       format('      data a(',i5,') /',i5,'/')
         if(a(i).ne.0)then
            write(k,59)i,i,i,hi(i),lo(i),xi(i)
59          format('      data hi(',i5,'),lo(',i5,'),xi(',i5,') /',
     $          i5,',',i5,',',1pe15.6,'/')
         end if
         do 4 j=1,vc
            write(k,60)j,i,c(j,i)
60          format('      data c(',i3,',',i5,') /',i5,'/')
    4    continue
    3 continue
      do 5 i=1,nv
         write(k,61)i,vval(0,i)
61       format('      data vval(0,',i5,') /',1pe15.6,'/')
         do 6 j=1,d
            write(k,62)i,j,v(i,j)
62          format('      data v(',i5,',',i2,') /',1pe15.6,'/')
            write(k,63)j,i,vval(j,i)
63          format('      data vval(',i2,',',i5,') /',1pe15.6,'/')
    6    continue
    5 continue
      write(k,*)'      loeval=ehg128(z,d,nc,vc,a,xi,lo,hi,c,v,nv,vval)'
      write(k,*)'      end'
      return
      end
      subroutine lofort(iunit,iv,liv,lv,wv)
      integer execnt,iunit
      integer iv(*)
      DOUBLE PRECISION wv(*)
      external ehg170
      save execnt
      data execnt /0/
      execnt=execnt+1
      call ehg170(iunit,iv(2),iv(4),iv(6),iv(14),iv(5),iv(17),iv(iv(7)),
     +iv(iv(8)),iv(iv(9)),iv(iv(10)),wv(iv(11)),wv(iv(13)),wv(iv(12)))
      return
      end
