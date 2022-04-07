Function ComputePmax,p,dth


x=indgen(90./dth+1)*dth
tmptabmaxp=fltarr(360./dth)*dth
for k=0,360./dth-1 do begin
y=gaussfit(x,p(*,k))
p2=p(*,k)-y

tmptabmaxp(k)=max(p2)
endfor

return,tmptabmaxp
end
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
Function ComputeInt,p,tabimin,dth,totalint  ; compute integral of sides lobes after removing first lobe
; -------------------------------------------------------------------
dthdph=(dth*!dtor)^2
grad=0.

x=indgen(90./dth+1)*dth
flagnonull=0
if (1) then begin
for k=0,360./dth-1 do begin          ; loop on azimut
   df=deriv(x,p(*,k))           ; compute 1D derivate
   dfprod=df*shift(df,1)        ; product of consecutive term, look for negative term (cross 0)
   w=where(dfprod(1:*) lt 0)
 
; TEST ZONE
wmax=where(abs(df) eq max(abs(df)))
df1=df(wmax(0):*)
df2=shift(df1,-1)
df3=df2-df1
df4=shift(df3,1)

wplat=where((df3*df4)(1:*) lt 0)
wfin=wmax(0)+wplat(0)+1
; END TESTZONE


  if wplat(0) ne -1 and w(0) ne -1 then begin     ; if found store location
      tabimin(k)=min([w(0),wfin(0)])
      flagnonull=0
   endif else if wplat(0) eq -1 and w(0) ne -1 then begin
      tabimin(k)=w(0)
   endif else if w(0) eq -1 and wplat(0) ne -1 then begin
      tabimin(k)=wfin(0)
   endif else if wplat(0) eq -1 and w(0) eq -1 then begin
        tabimin(k)=long(90/dth)
     endif

endfor
endif
if (0) then begin
for k=0,360./dth-1 do begin          ; loop on azimut
   df=deriv(x,p(*,k))           ; compute 1D derivate
   df2=deriv(x,df)
   dfprod=df*shift(df,1)        ; product of consecutive term, look for negative term (cross 0)
   dfprod2=df2*shift(df2,1)
   w=where(dfprod(1:*) lt 0)
   w2=where(dfprod2(1:*) lt 0)
   
  if w(0) ne -1 and w2(0) ne -1 then begin     ; if found store location
      tabimin(k)=min([w(0),w2(0)])
      flagnonull=0
   endif else if w2(0) eq -1 then begin
      tabimin(k)=w(0)
   endif else if w(0) eq -1 then begin
      tabimin(k)=w2(0)
   endif else if w(0) eq -1 and w2(0) eq -1 then begin
        tabimin(k)=long(90/dth)
     endif

endfor
endif

if (0) then begin
for k=0,360./dth-1 do begin          ; loop on azimut
   df=deriv(x,p(*,k))           ; compute 1D derivate
   df2=deriv(x,df)
   dfprod=df*shift(df,1)        ; product of consecutive term, look for negative term (cross 0)
   w=where(dfprod(1:*) lt 0)
   if w(0) ne -1 then begin     ; if found store location
      tabimin(k)=w(0)
      flagnonull=0
   endif else begin             ; if not found, look at 2nd derivate null
      
      dfprod2=df2*shift(df2,1)
      w2=where(dfprod2(1:*) lt 0)
    if w2(0) ne -1 then begin ; if found store location
        tabimin(k)=w2(0)
      endif else begin
        tabimin(k)=long(90/dth)
      endelse
   endelse

endfor
endif



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
totalint=total(p*tabsin*dthdph)

return,tmpInt
end

; -------------------------------------------------------------------
Function ComputeInt2,p,tabimin,dth,totalint  ; compute integral of sides lobes after removing first lobe
; -------------------------------------------------------------------
dthdph=(dth*!dtor)^2
grad=0.

x=indgen(90./dth+1)*dth
flagnonull=0
for k=0,360./dth-1 do begin          ; loop on azimut
   df=deriv(x,p(*,k))           ; compute 1D derivate
   dfprod=df*shift(df,1)        ; product of consecutive term, look for negative term (cross 0)
   w=where(dfprod(1:*) lt 0)
   if w(0) ne -1 then begin     ; if found store location
      tabimin(k)=w(0)
      flagnonull=0
;   endif ;else begin             ; if not found, look at 2nd derivate null
    ;  df2=deriv(x,df)
  ;    dfprod2=df2*shift(df2,1)
  ;    w2=where(dfprod2(1:*) lt 0)
 ;     if w2(0) ne -1 then begin ; if found store location
 ;        tabimin(k)=w2(0)
      endif else begin
         tabimin(k)=long(90/dth)
   ;   endelse
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
totalint=total(p*tabsin*dthdph)
return,tmpInt
end


; -------------------------------------------------------------------
  pro ARRAY_BEAM_SA_repuls2, numtry,n, a,k, aout,tabout,ENV=ENV, XCOS=XCOS, PLOT=PLOT, DB=DB
; -------------------------------------------------------------------
; MOD. J. Girard - 07/2012
;
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


    dth=2.
    aa=fltarr(n,3)
    aa0=fltarr(n,3)             ; save initial pattern
    aa=transpose(a)
    aa0=aa
    tmpa=fltarr(n,3)            ; New position array for trial
    insideplot=1

                                ; Stop parameters
    epsilon=1.e-7               ; epsilon of stopping criterion
    muststop=0                  ; Flag of stopping
    MAXITER=50000               ; Maximum iteration number

    distmin=2*sqrt(1./!pi/k)     ; minimum distance between antennas
    seed=long(systime(/seconds)) ; for antenna collision
    
    p=fltarr(90./dth+1,360./dth+1) ; power diagram

    tabDxn=fltarr(n,3)             ; Array of small moves
    tabimin=intarr(360./dth+1)     ; Array of first nulls

    tabInt=fltarr(1)            ; Array of Integral value (cnt,SLint,TotInt)
    tabtotint=fltarr(1)
    tabsll=fltarr(1)
    tabdist=fltarr(n,n)         ; Array of distances between antennas
    tabtemp=fltarr(1)           ; Array of temperature
    imin=0                      ; zenith angle of first null

    cnt=long(0)                       ; Loop counter
    cnt2=long(0)
    maxp=0
    maxp0=0
    int0=0
    int=0 
meanint=42.
meanmaxp=42.

    ; SA PARAMETERS
    Temp0=40000000              ; System initial temperature (HIGH)
    Temp=Temp0 
    Tfinal=0                    ; arbitrary final temperature
    Tstep=0.99                 ; Step in Cooling Schedule
    meanscore=0                 ; mean of cost function to compute initial temp
    kb=1; 1e-3                        ; Boltzmann constant analog
    G0=0.01 ; 0.0025                     ; Gain (in lambda unit)
    G=G0
    int=1e10                    ; Arbitrary high SL integral value


; Switch (0/1)
    noreheat=0
    aleacount=0                   ; Control counter to allow reheating
   
    name='saveSA'
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
    
    Ni=long(90./dth)+1          ; Define Power diagram float arrays
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
    

    p=ComputeP(aa,Ni,Nj,sc,ss,cc) ; Compute first PD
    tabint(0,0)=0
    
    p=p*e                           ; Apply enveloppe
    
    tabInt(0)=ComputeInt(p,tabimin,dth,ttint)/ttint ; Compute and store first Int
    int=tabint(0)
    print,"int0=",int
    tabtotint(0)=ttint
    tabsll(0)=0
    tabtemp(0)=Temp             ; Store first temperature
    p0=10*alog10(p)>0           ; Save initial PD
    cntreheat=0
; Here comes the simulated annealing code
    
    WHILE (Temp gt Tfinal and Muststop eq 0 and cntreheat lt 30) do begin
       cnt=cnt+1
       
       EQUIL=0                  ; Flag for thermodynamic equilibrium
       TDcount=0                
       tabtempint=0
       while EQUIL ne 30 do begin ; while TD equilibrium is not reached
          TDcount=TDcount+1
          aleacount=0                       ; Loop over random tries
          tabDxn(*,0)=randomu(seed,n)*G-G/2 ; random shifts centered on 0
          tabDxn(*,1)=randomu(seed,n)*G-G/2
          ArrdistAnt=fltarr(n,n)
          
                                ; Check if antennas are too close to each other
          REPEAT begin 
             
             tmpa(*,0)=aa(*,0)+tabDxn(*,0)
             tmpa(*,1)=aa(*,1)+tabDxn(*,1)
             
             for q1=0,n-1 do for q2=0,n-1 do begin
                Arrdistant(q1,q2)=sqrt((tmpa(q1,0)-tmpa(q2,0))^2+(tmpa(q1,1)-tmpa(q2,1))^2)
             endfor
             BadAntLoc=where(ArrDistAnt LT distmin and ArrDistAnt gt 0)
             
             if total(badantloc) ne -1. then begin
                for z=0,n_elements(BadAntLoc)-1 do begin
                   index=array_indices(arrdistant,badantloc)
                   ;angle=randomu(seed)*2*!pi
                   angle=atan(aa(index(0),1)-aa(index(1),1),aa(index(0),0)-aa(index(1),0))
                   tabdxn(index(0),0)=tabdxn(index(0),0)+distmin*1/50*cos(angle) ; move one in x
                   tabdxn(index(1),0)= tabdxn(index(1),0)-distmin*1/50*cos(angle)
                   tabdxn(index(0),1)=tabdxn(index(0),1)+distmin*1/50*sin(angle)
                   tabdxn(index(1),1)=tabdxn(index(1),1)-distmin*1/50*sin(angle) ; move in y
                                ; rq: due to matrix symetry => dmin/2   <<== it sucks !to improve!
                endfor
             endif
          ENDREP until (total(ArrdistAnt LT distmin and ArrDistAnt gt 0) eq 0.)
          

; If system is stuck (cost function not minimized since the last n tests)
          if TDcount gt 2000 and noreheat eq 0 then begin
                temp=temp*(1/Tstep)^10   ; Raise current temperature (move n steps backward)
               ; ; G=G*0.9                ; modify gain
               
                print,'REHEATING'
                TDcount=0
                noreheat=1               ; disable reheating
                cntreheat=cntreheat+1
             endif

; Compute new power diagram and store integral value of SL from first null
          newp=ComputeP(aa+tabDxn,Ni,Nj,sc,ss,cc)
          newp=newp*e           ; apply antenna enveloppe
          
          tabimin=intarr(360./dth)
          Int=ComputeInt(newp,tabimin,dth,ttint) ; Compute the value of the cost function
          ;print,int,ttint
          int=int/ttint
          ;print,int
        
; look for Worst Side lobe just for curiosity (nothing to do with
; Kogan)
         ; tabmaxp=computepmax(newp,dth)
         ; tabmaxp=10*alog10(tabmaxp/max(newp))<0
          newp=(10.*alog10(newp/max(newp)))<0     ; to normalized BP in dB
          tabmaxp=fltarr(n_elements(tabimin)-1)
          for kimin=0,n_elements(tabimin)-2 do tabmaxp(kimin)=max(newp(min(tabimin):*,kimin)) ; fill sidelobe levels for each phi
          locmax=(where(tabmaxp eq max(tabmaxp)))(0) ; give phi of worst side lobe
         ; maxP=max(newp(min(tabimin):*,*)) ; force to choose local or absolute phi of first null
 
          maxP=max(tabmaxp);newp(tabimin(locmax):*,locmax)) ; force to choose local or absolute phi of first null

          if maxp lt maxprecord then begin
             maxprecord=maxp
             aarecord=aa
          endif
        
          w=where(tabimin gt 0)
          if total(w) ne -1 then begin
          imin=min(tabimin(w))
       endif else imin=min(tabimin)
          locimin=(where(tabimin eq imin))(0)
 
          Score=int-tabint(cnt-1) ; Compare new "energy" value with previous value
          ;print,"score=",score,int
         ; delta=score/int         ; Relative variation of CF
;print,"Boltzmann threshold=",exp(-Score/(kb*Temp))*100<100.,"%"
        
          if Score lt 0 then begin ; if a state of lower energy has been found
             aa=aa+tabDxn          ; keep it!
             equil=equil+1         ; get closer to equilibrium
             noreheat=0            ; enable reheating
            ; print,'normal'
            ; if (abs(delta) lt 0.1) then EQUIL=EQUIL+1 else EQUIL=0 ; if system start to reach equilibrium
      ;  print,'DCF=',score,' kbT=',kb*Temp,' ratio=',-score/(kb*temp),' ==> Steepest ..'
          endif else if (randomu(seed) lt exp(-Score/(kb*Temp))) then begin ; else try to go through the potential barrier
             aa=aa+tabDxn  
            ; print,'chanceux!'                   ; if succeeded then keep new config
;print,'DCF=',score,' kbT=',kb*Temp,' ratio=',-score/(kb*temp)
            ;if (abs(delta) lt 0.1) then EQUIL=EQUIL+1 else EQUIL=0
             equil=equil+1         ; get closer to equilibrium
             noreheat=0            ; enable reheating
          endif

cnt2=cnt2+TDcount
endwhile



tabtotint=[tabtotint,ttint]
      tabInt=[tabint,Int] ;store integral value
      tabtemp=[tabtemp,temp]
      tabsll=[tabsll,maxp]

      
      
     
         if (cnt eq MAXITER) then begin
         Muststop = 1
            print,"System is now cristallized!"
         endif else begin 
            Muststop=0
            print,"System is cooling..."
         endelse
     ; endif

       Temp=Temp*Tstep          ; Let's cool the system a bit, but slowly!
       save,file=name,aa0,p0,aa,p,tabint,tabtotint,tabtemp ; save current state, just in case...
       if (cnt MOD 100 eq 0) then print,cnt       ; print cnt every 100 steps

; Measure system temperature
       if (cnt lt 20 and cnt gt 2) then begin
print,"score=",score
meanscore=[meanscore,[score]] 
endif

       if cnt eq 20 then begin
          Temp=mean(abs(meanscore))
       endif

       if cnt eq 1 then begin
          tabmeanmaxp=maxp
          tabmeanint=int
       endif

       if cnt lt 10 then begin 
          maxp0=maxp0+maxp/10
          int0=int0+int/10
       endif

       if cnt le 100  and cnt gt 1 then begin 
          tabmeanmaxp=[tabmeanmaxp,[maxp]]
          tabmeanint=[tabmeanint,[int]]
          print,int,maxp,ttint
       endif

       if cnt gt 100 then begin 
          tabmeanmaxp=[tabmeanmaxp(1:*),[maxp]]
          tabmeanint=[tabmeanint(1:*),[int]]
          meanmaxp=total(tabmeanmaxp)/n_elements(tabmeanmaxp)
          meanint=total(tabmeanint)/n_elements(tabmeanint)
          sigman=sigma(tabmeanint)/meanint
          
          if abs(sigman) lt epsilon then begin
             muststop=1
             print,"under epsilon ==> stopping at next iter"
          endif
   
          print,'cnt=',cnt,'cnt2=',cnt2,'G=',G,'int=',int,' sll=',maxp,' sigma',sigman,' temp',temp,' kb*temp',kb*temp,' score',score,"Boltzmann threshold=",exp(-Score/(kb*Temp))*100<100.,"%"
       endif

; PLOT
  if (InsidePlot eq 1 and (cnt MOD 10 eq 0 or cnt eq 1)) then begin ; Artist mode (poetry and paintings)
        
; unreadable plot magic spells
          ; Plot antenna pattern with effective area       
          usersym,sin(findgen(17)*!pi/8),cos(findgen(17)*!pi/8),/fill
          plot,aa(*,0),aa(*,1),xra=[xmin,xmax],yra=[ymin,ymax],/xsty,/ysty,psym=8,title='SA Antenna pattern Dmin=0.39',xtit='Wavelengths',ytit='Wavelengths',/isotropic,color=0,background=250
          for i=0,n-1 do oplot,aa(i,0)+cos(findgen(181)*2*!dtor)*sqrt(1./!pi/k),aa(i,1)+sin(findgen(181)*2*!dtor)*sqrt(1./!pi/k),color=0

          ; Plot power diagram
          
          plot,the,newp(*,0),xr=[0,90],yr=[-50,0],/xs,/ys,color=0,background=250,/nodata,title='Power Diagram',xtit='ZA (deg.)',ytit="Pow (dB)"

          for j=1,360./dth do oplot,the,newp(*,j),color=0
          oplot,the,newp(*,imin),color=100,thick=2

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
          xyouts,72,-15,"Az = "+string(imin),color=0,charsize=1.5
      
     
    endif


       
    ENDWHILE


Mmaxp=meanmaxp
Mint=meanint
print,cnt
tabout=[imin,maxp,maxp0,Mmaxp,int,int0,Mint,maxprecord,cnt]
aout=transpose(aarecord)

save,file="savesa-"+strcompress(string(n),/remove_all)+'-'+strcompress(string(numtry),/remove_all),aa0,aa,aarecord,p0,newp,tabint,tabtemp,tabsll,tabtotint
aout=transpose(aarecord)

    return
 end
  
  
