; -------------------------------------------------------------------
Function ComputeP,aa,Ni,Nj,sc,ss,cc ; compute Power diagram given antenna pattern aa ( n * 3 )
; -------------------------------------------------------------------
;common champ, field  ;global, needed in ComputeInt
field=fltarr(Ni,Nj)
power=fltarr(Ni,Nj)
eiphi=complex(0.,0.)
; Compute P(i,j)
;for i=0,45 do for j=0,180 do begin
;   th=(i*2)*!dtor               ; 0,90 deg, step 2 deg
;   ph=(j*2)*!dtor               ; 0,360 deg, step 2 deg
;   u=[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]
;   phi=2*!pi*(aa#u)
;   eiphi=total(complex(cos(phi),sin(phi)))
;   power(i,j)=float(eiphi*conj(eiphi))
;   field(i,j)=real_part(eiphi)
;endfor

n=(size(aa))(1)
for i=0,N-1 do begin
phi=2*!pi*(aa(i,0)*sc+aa(i,1)*ss+aa(i,2)*cc)
eiphi=eiphi+complex(cos(phi),sin(phi))
endfor


;field=real_part(eiphi)
power=float(eiphi*conj(eiphi))

return,power
end

; -------------------------------------------------------------------
Function ComputeInt,p,tabimin,dth  ; compute integral of sides lobes after removing first lobe
; -------------------------------------------------------------------
dthdph=(dth*!dtor)^2
grad=0.

x=indgen(90./dth+1)*dth
flagnonull=0
for k=0,360./dth do begin          ; loop on azimut
   df=deriv(x,p(*,k))           ; compute 1D derivate
   dfprod=df*shift(df,1)        ; product of consecutive term, look for negative term (cross 0)
   w=where(dfprod(1:*) lt 0)
   if w(0) ne -1 then begin     ; if found store location
      tabimin(k)=w(0)
      flagnonull=0
 ;  endif else begin             ; if not found, look at 2nd derivate null
 ;     df2=deriv(x,df)
 ;     dfprod2=df2*shift(df2,1)
  ;    w2=where(dfprod2(1:*) lt 0)
   ;   if w2(0) ne -1 then begin ; if found store location
   ;      tabimin(k)=w2(0)
      endif else begin
         tabimin(k)=long(90./dth)
    ;  endelse
   endelse
endfor

tmpInt=0.0

w=where(tabimin gt 0)
if total(w) ne -1 then begin
   imin=min(tabimin(w))     ; first encoutered null (in zenith angle direction) in first null array
endif else imin=min(tabimin)

; Compute omega_SL integral (after 1st null)
; integration over theta,phi
tabsin=sin(!dtor*dth*findgen(90./dth+1))
tabsin=rebin(reform(tabsin,90./dth+1,1),90./dth+1,360./dth+1)
;tmpInt=tmpInt+total(p(min(tabimin(w)):45,*)*tabsin(min(tabimin(w)):45,*)*dthdph)
tmpInt=tmpInt+total(p(imin:long(90./dth),*)*tabsin(imin:long(90./dth),*)*dthdph) ; integrate from first null
;tmpInt=tmpInt+total(p(min(tabimin):long(90./dth),*)*tabsin(min(tabimin):lon(90./dth),*)*dthdph) ; integrate from "first" first null
;stop
return,tmpInt
end

; -------------------------------------------------------------------
  pro ARRAY_BEAM_KOGANFIN, numtry,n, a, k,aout,tabout,ENV=ENV, XCOS=XCOS, PLOT=PLOT, DB=DB
; -------------------------------------------------------------------
;
; MOD. J. Girard - 07/2012
; [IN]     numtry = number of antennas
; [IN]     n = number of antennas
; [IN]     a(3,n) = (x,y) coordinates of antennas (normalized to wavelength)
; [IN]     lambda^2 / k = effective area of elementary antenna
;          (k=8 -> Dipole, k=3 -> DAM, ...)
; [OUT]    p = power diagram
; [OUT]    aa = current antenna position
; [KEY]    ENV = envelope (beam of elementary antenna) (46) or (46,181)
; [KEY]    XCOS = envelope in cos(theta)^XCOS (2 = Dipole, 1 ~ DAM, ...)
; [KEY]    PLOT = CARTESIAN, POLAR, MAP  display
; [KEY]    DB = display in DB (default = linear)
  

  

  print,"Number of Antennas=",n," Try number=",numtry
    dth=2.
    aa=fltarr(n,3)
    aa0=fltarr(n,3)       
    aa=transpose(a)
    aa0=aa

    epsilon=1e-7                ; epsilon of stopping criterion
    muststop=0
    MAXITER=50000
    k=8
    distmin=2*sqrt(1./!pi/k)    ; minimum distance between antennas
insideplot=0

if insideplot eq 1 then begin
 window,0,xs=900,ys=450      
   device,decomposed=0
   TVLCT,[[0],[0],[0]],0       ;black
   TVLCT,[[255],[0],[0]],100   ;R
   TVLCT,[[0],[255],[0]],150   ;V
   TVLCT,[[0],[0],[255]],200   ;B
   TVLCT,[[255],[255],[255]],250 ;Wh
   erase
   !p.multi=[0,2,1,0,0]
   !p.charsize=2

xmin=-2.
xmax=2.
ymin=xmin
ymax=xmax
endif

    p=fltarr(90./dth+1,360./dth+1)           ; power diagram
   
    tabDxn=fltarr(n,3)
    tabimin=intarr(360./dth+1)

    tabInt=fltarr(1)
    tabsll=fltarr(1)
    tabdist=fltarr(n,n)         ; Array of distances between antennas
    tabtemp=fltarr(1)           ; Array of temperature
    imin=0

    cnt=long(0)
    maxp=0
    maxp0=0
    int0=0
    int=0 


    ; KOGAN PARAMETERS
    G=0.010
    newWSL=0
    iloc=0
    jloc=0

    R=2 ;?


    name='saveKO'

; Enveloppe
    the=findgen(90./dth+1)*2    ; Zenith angle
    cth=abs(cos(the*!dtor))
    sth=abs(sin(the*!dtor))
    e=fltarr(90./dth+1,360./dth+1)+1. ; enveloppe of antenna

     if keyword_set(ENV) then begin
       e=size(ENV)
       if e(0) eq 2 then e=ENV else e=rebin(reform(ENV,90./dth+1,1),90./dth+1,360./dth+1)
    endif else if keyword_set(XCOS) then e=rebin(reform(cth^XCOS,90./dth+1,1),90./dth+1,360./dth+1)


aarecord=aa
maxprecord=maxp 

Ni=long(90./dth)+1              ; Define Power diagram float array
Nj=long(360./dth)+1
i=1.*indgen(Ni)*dth
j=1.*indgen(Nj)*dth
th=rebin(i,Ni,Nj)
ph=transpose(rebin(j,Nj,Ni))
th=th*!dtor
ph=ph*!dtor

sc=sin(th)*cos(ph)
ss=sin(th)*sin(ph)
cc=cos(th)

p=ComputeP(aa,Ni,Nj,sc,ss,cc)   ; Compute first PD
tabint(0,0)=0
p=p*e                       ; Apply enveloppe
tabInt(0)=ComputeInt(p,tabimin,dth) ; Compute and store first Int
tabsll(0)=0
p0=10*alog10(p)>0                        ; Save initial PD

; Here comes KOGAN code

REPEAT begin
   newp=ComputeP(aa,Ni,Nj,sc,ss,cc)
   newp=newp*e                  ; applying the antenna response pattern
   Int=ComputeInt(newp,tabimin,dth) 
                                         
   newp=(10.*alog10(newp/max(newp)))<0     ; to normalized BP in dB
  ; tvscl,newp
   
   ;wait,0.3
; select the worst residual lobe profile (max of residual integral)
   tabmaxp=fltarr(n_elements(tabimin)-1)
   for kimin=0,n_elements(tabimin)-2 do tabmaxp(kimin)=max(newp(min(tabimin):*,kimin))
   locmax=(where(tabmaxp eq max(tabmaxp)))(0) ;=> give phi
   maxP=max(newp(min(tabimin):*,locmax))
   LocWSL=(where(newp eq maxP))(0)
    ;print,cnt,maxp,LocWSL
   ; stop
   if maxp lt maxprecord then begin
      maxprecord=maxp
      aarecord=aa
   endif

      tabInt=[tabint,Int]
      tabsll=[tabsll,maxP] 

      index=array_indices(newp,LocWSL) ; location of WSL in newp
      iLoc=index(0)
      jloc=index(1)          
      Wph=(jloc*dth)*!dtor       
      Wth=(iloc*dth)*!dtor
      
      w=where(tabimin gt 0)
          if total(w) ne -1 then begin
          imin=min(tabimin(w))
       endif else imin=min(tabimin)
          locimin=(where(tabimin eq imin))(0)

     ; locimin=(where(tabimin eq min(tabimin)))(0)
;wait,0.3

; compute Dxn for each antenna and return x and y projection in tabDxn
; tabr=sqrt(aa(*,0)^2+aa(*,1)^2)
      for k1=0,n-1 do begin 
         dA=0.
         dA=total(sin(2*!pi*(sin(wth)*cos(wph)*(aa(*,0)-aa(k1,0))+sin(wth)*sin(wph)*(aa(*,1)-aa(k1,1)))))
         tabDxn(k1,0)=-G*dA*sin(wth)*cos(wph)
         tabDxn(k1,1)=-G*dA*sin(wth)*sin(wph)
      endfor

      aa(*,0)=aa(*,0)+tabdxn(*,0)
      aa(*,1)=aa(*,1)+tabdxn(*,1)

; Control if distrution meets constraints (distance(k,k') >= distmin m and
; within a circle of radius R)
      
      ArrdistAnt=fltarr(n,n)

; Check if antennas are too close to each other
      REPEAT begin 
         
         for q1=0,n-1 do for q2=0,n-1 do begin
            Arrdistant(q1,q2)=sqrt((aa(q1,0)-aa(q2,0))^2+(aa(q1,1)-aa(q2,1))^2)
         endfor
         BadAntLoc=where(ArrDistAnt LT distmin and ArrDistAnt gt 0)
         
         if total(badantloc) ne -1. then begin
            for z=0,n_elements(BadAntLoc)-1 do begin
              ; print,"cor!"
               index=array_indices(arrdistant,badantloc)
               angle=atan(aa(index(0),1)-aa(index(1),1),aa(index(0),0)-aa(index(1),0))
               aa(index(0),0)=aa(index(0),0)+distmin*1./50*cos(angle) ; move one in x
               aa(index(1),0)=aa(index(1),0)-distmin*1./50*cos(angle)
               aa(index(0),1)=aa(index(0),1)+distmin*1./50*sin(angle)
               aa(index(1),1)=aa(index(1),1)-distmin*1./50*sin(angle) ; move in y
                                ; rq: due to matrix symetry => dmin/2   <<== it sucks !to improve!
            endfor
         endif
         
      ENDREP until (total(ArrdistAnt LT distmin and ArrDistAnt gt 0) eq 0.)
      
     
      if cnt eq 0 then begin
         tabmeanmaxp=maxp
         tabmeanint=int
         tabsigmamean=1
      endif

      if cnt lt 10 then begin 
         maxp0=maxp0+maxp/10
         int0=int0+int/10
      endif

      if cnt le 100  and cnt ge 1 then begin 
         tabMeanmaxp=[tabmeanmaxp,[maxp]]
         tabMeanint=[tabmeanint,[int]]
      endif
      
      if cnt gt 100 then begin 
         tabmeanmaxp=[tabmeanmaxp(1:*),[maxp]]
         tabmeanint=[tabmeanint(1:*),[int]]
         
         meanmaxp=mean(tabmeanmaxp) ;total(tabmeanmaxp)/n_elements(tabmeanmaxp)
         meanint=mean(tabmeanint) ;total(tabmeanint)/n_elements(tabmeanint)
         sigman=sigma(tabmeanmaxp)/meanmaxp
         if cnt gt 100 and cnt le 200 then begin
            tabsigmamean=[tabsigmamean,[sigman]]
            endif else if cnt gt 200 then begin
            tabsigmamean=[tabsigmamean(1:*),[sigman]]
            endif

         sigsig=sigma(tabsigmamean)/mean(tabsigmamean)
         if (cnt MOD 1000 eq 0) then print,sigsig,epsilon
         if abs(sigsig) lt epsilon then begin
            muststop=1
            print,cnt,abs(sigsig),epsilon
            
            print,"under epsilon ==> stopping at next iter"
         endif
       ;  print,'G=',G,'int=',int,' sll=',maxp,' sigma',sigman
      endif

if (cnt MOD 1000 eq 0) then print,cnt
cnt=cnt+1

; PLOT
  if (InsidePlot eq 1 and (cnt MOD 10 eq 0 or cnt eq 1)) then begin ; Artist mode (poetry and paintings)
        
; unreadable plot magic spells
          ; Plot antenna pattern with effective area       
          usersym,sin(findgen(17)*!pi/8),cos(findgen(17)*!pi/8),/fill
          plot,aa(*,0),aa(*,1),xra=[xmin,xmax],yra=[ymin,ymax],/xsty,/ysty,psym=8,title='SA Antenna pattern Dmin=0.39',xtit='Wavelengths',ytit='Wavelengths',/isotropic,color=0,background=250
          for i=0,n-1 do oplot,aa(i,0)+cos(findgen(181)*2*!dtor)*sqrt(1./!pi/k),aa(i,1)+sin(findgen(181)*2*!dtor)*sqrt(1./!pi/k),color=0

          ; Plot power diagram
          
          plot,the,newp(*,0),xr=[0,90],yr=[-50,0],/xs,/ys,color=0,background=250,/nodata,title='Power Diagram',xtit='ZA (deg.)',ytit="Pow (dB)"

          for j=1,180 do oplot,the,newp(*,j),color=0
          oplot,the,newp(*,locmax),color=100,thick=2

          oplot,[imin*dth,imin*dth],[-50,50],color=150,thick=2
          oplot,[0,90],[maxP,maxP],color=0

          ; labels
          xyouts,60,-7.5,strcompress("Loop "+string(cnt)),color=0,charsize=1.5
          xyouts,60,-10,strcompress("Int="+string(int)),color=0,charsize=1.5
          xyouts,20,-5,strcompress('2nd Max='+string(maxP)+" dB"),color=0,charsize=1.5
        ;  if countplotreheat ne 0 then begin       
        ;     xyouts,[20,20],[maxP,maxP],strcompress("REHEATING!"+string(countplotreheat)),color=100,charsize=3
        ;     countplotreheat=countplotreheat-1
        ;  endif
          
          ; legend
          oplot,[60,70],[-15,-15],color=100,thick=3
          xyouts,72,-15,"Az = 0",color=0,charsize=1.5
      
     
    endif



ENDREP until (cnt eq MAXITER or muststop eq 1)



save,file="saveko-"+strcompress(string(n),/remove_all)+'-'+strcompress(string(numtry),/remove_all),aa0,aa,aarecord,p0,newp,tabint,tabsll

Mmaxp=meanmaxp
Mint=meanint
print,cnt
tabout=[locimin,maxp,maxp0,Mmaxp,int,int0,Mint,maxprecord,cnt]
aout=transpose(aarecord)

return
end
  
