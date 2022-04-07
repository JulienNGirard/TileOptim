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
for k=0,360./dth-1 do begin          ; loop on azimut
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
Function ComputeInt2,p,tabimin,dth,totalint  ; compute integral of sides lobes after removing first lobe
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

function computepowdiag,na,ntry,aa,dth,newp,stat,ENV=ENV
; Preparation Computing Pmax
dthdph=(dth*!dtor)^2
x=indgen(90./dth+1)*dth

; Enveloppe
    the=findgen(90./dth+1)*2    ; Zenith angle
    cth=abs(cos(the*!dtor))
    sth=abs(sin(the*!dtor))
    e=fltarr(90./dth+1,360./dth+1)+1. ; enveloppe of antenna
    
    if keyword_set(ENV) then begin
       print,"Enveloppe!"
       xcos=2
       e=rebin(reform(cth^XCOS,90./dth+1,1),90./dth+1,360./dth+1)
    endif


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
    
    
    newp=ComputeP(aa,Ni,Nj,sc,ss,cc)
    newp=newp*e                 ; apply antenna enveloppe

    tabimin=intarr(360./dth)
    int=ComputeInt(newp,tabimin,dth);,ttint) ; Compute and store first Int
    newp=10*alog10(newp/max(newp))

    stat=computeFWHM(newp,dth)

    
    tabmaxp=fltarr(n_elements(tabimin)-1)
   ; w=where(tabimin gt 0 and tabimin lt 45)
 
for kimin=0,n_elements(tabimin)-2 do begin 
locm=min(tabimin)+2<long(90./dth)
tabmaxp(kimin)=max(newp(locm:*,kimin)) ; fill sidelobe levels for each phi
endfor


    locmax=(where(tabmaxp eq max(tabmaxp)))(0)                                          ; give phi of worst side lobe
    maxP=max(tabmaxp);newp(tabimin(locmax):*,locmax)) ; force to choose local or absolute phi of first null
    
    w=where(tabimin gt 0); and tabimin lt 45)
    if total(w) ne -1 then begin
       imin=min(tabimin(w))
    endif else imin=min(tabimin)
    locimin=(where(tabimin eq imin))(0)
    
    ; Plotting
   plot,the,newp(*,0),xr=[0,90],yr=[-80,0],/xs,/ys,color=0,background=250,/nodata,xtit='ZA (deg.)',ytit="Pow (dB)";,title='Power Diagram'

   for j=1,360./dth do oplot,the,newp(*,j),color=0
   oplot,the,newp(*,imin),color=100,thick=2
   
   oplot,[imin*dth,imin*dth],[-80,80],color=150,thick=2
   oplot,[0,90],[maxP,maxP],color=0
   oplot,[stat(0)/2,stat(0)/2],[-100,0],color=0
                                ; labels
   ;xyouts,60,-7.5,strcompress("Loop "+string(cnt)),color=0,charsize=1.5
  ; xyouts,60,-10,strcompress("Int="+string(int)),color=0,charsize=1.5
   xyouts,20,-5,strcompress('2nd Max='+string(maxP)+" dB"),color=0,charsize=1.5
   xyouts,20,-10,strcompress('FWHM='+string(stat(0))+" deg"),color=0,charsize=1.5
                                ; legend
  ; oplot,[60,70],[-15,-15],color=100,thick=3
  ; xyouts,72,-15,"Az = "+string(imin),color=0,charsize=1.5

return,maxp

end
