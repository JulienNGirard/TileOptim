pro plotres
; usage:
Ndeb=5
Nfin=22
Ntry=10
dn=Nfin-Ndeb+1
liste=indgen(dn)+ndeb

namerec="REC-"


restore,file='BestKOGSLL.sav'
restore,file='BestSASLL.sav'
restore,file='BestSISLL.sav'

; Loading savesets
restore,file='MES-KOG.sav'
MsllKO=tabmesmaxp
MFWHMKO=tabmesFWHM
restore,file='MES-'+namerec+'KOG.sav'
MsllKOr=tabmesmaxp
MFWHMKOr=tabmesFWHM

restore,file='MES-SA.sav'
MsllSA=tabmesmaxp
MFWHMSA=tabmesFWHM
restore,file='MES-'+namerec+'SA.sav'
MsllSAr=tabmesmaxp
MFWHMSAr=tabmesFWHM

restore,file='MES-SI.sav'
MsllSI=tabmesmaxp
MFWHMSI=tabmesFWHM
restore,file='MES-'+namerec+'SI.sav'
MsllSIr=tabmesmaxp
MFWHMSIr=tabmesFWHM

; loading masks

maskKOGr=tabKOGmaskrec
maskSAr=tabSAmaskrec
maskSIr=tabSImaskrec

maskKOG=tabKOGmaskend
maskSA=tabSAmaskend
maskSI=tabSImaskend

MsllKO2=fltarr(dn)
MsllSA2=fltarr(dn)
MsllSI2=fltarr(dn)
MsllKO2r=fltarr(dn)
MsllSA2r=fltarr(dn)
MsllSI2r=fltarr(dn)

MsllKO3=fltarr(dn)
MsllSA3=fltarr(dn)
MsllSI3=fltarr(dn)

if (1) then begin

statko=fltarr(dn,ntry)
statsa=fltarr(dn,ntry)
statsi=fltarr(dn,ntry)

tabstatko=fltarr(dn,3)
tabstatsa=fltarr(dn,3)
tabstatsi=fltarr(dn,3)


statko_f=fltarr(dn,ntry,6)
statsa_f=fltarr(dn,ntry,6)
statsi_f=fltarr(dn,ntry,6)

tabstatko_FWHM=fltarr(dn,3)
tabstatsa_FWHM=fltarr(dn,3)
tabstatsi_FWHM=fltarr(dn,3)

tabstatko_AR=fltarr(dn,3)
tabstatsa_AR=fltarr(dn,3)
tabstatsi_AR=fltarr(dn,3)

listesolko=intarr(dn,ntry) ; 1 for end 2 for rec
listesolsa=intarr(dn,ntry)
listesolsi=intarr(dn,ntry)

listesolkosa=intarr(dn)
listesolminko=intarr(dn)
listesolminsa=intarr(dn)

for i=0,dn-1 do begin

for itry=0,ntry-1 do begin
;statko(i,itry)=msllKO(i,itry) gt msllKOr(i,itry)

   if msllKO(i,itry) gt -100. and msllKOr(i,itry) gt -100. then begin
      
      if msllKO(i,itry) lt msllKOr(i,itry) then begin
         statko(i,itry)=msllKO(i,itry)
         statko_f(i,itry,*)=MFWHMKO(i,itry,*)
         listesolko(i,itry)=1
      endif else begin
         statko(i,itry)=msllKOr(i,itry)
         statko_f(i,itry,*)=MFWHMKOr(i,itry,*)
         listesolko(i,itry)=2
      endelse
      
; statko(i,itry)=min([msllKO(i,itry),msllKOr(i,itry)])


   endif else if msllKO(i,itry) lt -100. then begin
      statko(i,itry)=msllKOr(i,itry)
      statko_f(i,itry,*)=MFWHMKOr(i,itry,*)
      listesolko(i,itry)=2
   endif else if msllKOr(i,itry) lt -100. then begin
      statko(i,itry)=msllKO(i,itry)
      statko_f(i,itry,*)=MFWHMKO(i,itry,*)
      listesolko(i,itry)=1
   endif

    if msllSA(i,itry) gt -100. and msllSAr(i,itry) gt -100. then begin

       if msllSA(i,itry) lt msllSAr(i,itry) then begin
          statsa(i,itry)=msllSA(i,itry)
          statsa_f(i,itry,*)=MFWHMSA(i,itry,*)
          listesolsa(i,itry)=1
       endif else begin
          statsa(i,itry)=msllSAr(i,itry)
          statsa_f(i,itry,*)=MFWHMSAr(i,itry,*)
          listesolsa(i,itry)=2
       endelse

   ;   statsa(i,itry)=min([msllSA(i,itry),msllSAr(i,itry)])
   endif else if msllSA(i,itry) lt -100. then begin
      statsa(i,itry)=msllsar(i,itry)
      statsa_f(i,itry,*)=MFWHMSAr(i,itry,*)
      listesolsa(i,itry)=2
   endif else if msllSAr(i,itry) lt -100. then begin
      statsa(i,itry)=msllsa(i,itry)
      statsa_f(i,itry,*)=MFWHMSA(i,itry,*)
      listesolsa(i,itry)=1
   endif

  if msllSI(i,itry) gt -100. and msllSIr(i,itry) gt -100. then begin
     
     if msllSA(i,itry) lt msllSAr(i,itry) then begin
        statsi(i,itry)=msllSI(i,itry)
        statsi_f(i,itry,*)=MFWHMSI(i,itry,*)
        listesolsi(i,itry)=1
     endif else begin
        statsi(i,itry)=msllSIr(i,itry)
        statsi_f(i,itry,*)=MFWHMSIr(i,itry,*)
        listesolsi(i,itry)=2
     endelse
     
; statsi(i,itry)=min([msllSI(i,itry),msllSIr(i,itry)])
   endif else if msllSI(i,itry) lt -100. then begin
      statsi(i,itry)=msllsir(i,itry)
      statsi_f(i,itry,*)=MFWHMSIr(i,itry,*)
      listesolsi(i,itry)=2
   endif else if msllSIr(i,itry) lt -100. then begin
      statsi(i,itry)=msllsi(i,itry)
      statsi_f(i,itry,*)=MFWHMSI(i,itry,*)
      listesolsi(i,itry)=1
   endif
endfor




wko=where(statko(i,*) lt -100.)
wsa=where(statsa(i,*) lt -100.)
wsi=where(statsi(i,*) lt -100.)

tabstatko(i,0)=mean(statko(i,*))
tabstatko(i,1)=sigma(statko(i,*))
tabstatko(i,2)=min(statko(i,*))
w=where(statko(i,*) eq min(statko(i,*)))

listesolminko(i)=w
;if total(wko) ne -1 then tabstatko(i,1)=0

tabstatsa(i,0)=mean(statsa(i,*))
tabstatsa(i,1)=sigma(statsa(i,*))
tabstatsa(i,2)=min(statsa(i,*))
w=where(statsa(i,*) eq min(statsa(i,*)))
listesolminsa(i)=w
;if total(wsa) ne -1 then tabstatsa(i,1)=0

tabstatsi(i,0)=mean(statsi(i,*))
tabstatsi(i,1)=sigma(statsi(i,*))
tabstatsi(i,2)=min(statsi(i,*))
;if total(wsi) ne -1 then tabstatsi(i,1)=0


tabstatko_FWHM(i,0)=mean(statko_f(i,*,0))
tabstatko_FWHM(i,1)=sigma(statko_f(i,*,0))
tabstatko_FWHM(i,2)=max(statko_f(i,*,0))

tabstatko_AR(i,0)=mean(statko_f(i,*,5))
tabstatko_AR(i,1)=sigma(statko_f(i,*,5))
tabstatko_AR(i,2)=min(statko_f(i,*,5))

tabstatsa_FWHM(i,0)=mean(statsa_f(i,*,0))
tabstatsa_FWHM(i,1)=sigma(statsa_f(i,*,0))
tabstatsa_FWHM(i,2)=max(statsa_f(i,*,0))

tabstatsa_AR(i,0)=mean(statsa_f(i,*,5))
tabstatsa_AR(i,1)=sigma(statsa_f(i,*,5))
tabstatsa_AR(i,2)=min(statsa_f(i,*,5))

tabstatsi_FWHM(i,0)=mean(statsa_f(i,*,0))
tabstatsi_FWHM(i,1)=sigma(statsa_f(i,*,0))
tabstatsi_FWHM(i,2)=max(statsa_f(i,*,0))

tabstatsi_AR(i,0)=mean(statsa_f(i,*,5))
tabstatsi_AR(i,1)=sigma(statsa_f(i,*,5))
tabstatsi_AR(i,2)=min(statsa_f(i,*,5))

endfor
endif

if (0) then begin
for i=0,dn-1 do begin
MsllKO2(i)=min(MsllKO(i,where(maskKOG(*,i) eq 1)))
MsllSA2(i)=min(MsllSA(i,where(maskSA(*,i) eq 1)))
MsllSI2(i)=min(MsllSI(i,where(maskSI(*,i) eq 1)))

MsllKO2r(i)=min(MsllKOr(i,where(maskKOGr(*,i) eq 1)))
MsllSA2r(i)=min(MsllSAr(i,where(maskSAr(*,i) eq 1)))
MsllSI2r(i)=min(MsllSIr(i,where(maskSIr(*,i) eq 1)))

MsllKO3(i)=min([MsllKO2r(i),MsllKO2(i)])
MsllSA3(i)=min([MsllSA2r(i),MsllSA2(i)])
MsllSI3(i)=min([MsllSI2r(i),MsllSI2(i)])
endfor
endif



set_plot,'ps'
device, yoff=1., ysize=15., xoff=0., xsize=20, /portrait
device,/color,bits_per_pixel=8
device,file='3-ALL-'+'MIN-'+'5-22.ps'
loadct,39
!p.multi=0
!p.font=0

; usersym
r=0.7
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=4 ;,/fill

; Plot initial distribution SLL
plot,listeKOG,tabKOGSLL(0,*),title='Lowest SLL = f(Nant)',xtit='Number of antennas',ytit='SLL (dB)',/nodata,yr=[-70,max(tabSASLL(0,*))],charsize=2,thick=4,xr=[4,23],/xsty
oploterror,listeKOG,tabKOGSLL(0,*),tabKOGSLL(1,*),psym=8,symsize=1,thick=4

; Fit a curve on it
;yfit=linfit(10*alog10(listeKOG),reform(tabKOGSLL(0,*)))
yfit=linfit(listeKOG,reform(tabKOGSLL(0,*)))

print,'Fit coefficient',' slope=',yfit(1),' y0=',yfit(0)
stop
;oplot,listeKOG,yfit(1)*10*alog10(listeKOG)+yfit(0),thick=4,line=1
oplot,listeKOG,yfit(1)*(listeKOG)+yfit(0),thick=4,line=1
; usersym
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=1,/fill

if (0) then begin
oplot,listeKOG,MsllKO3,psym=8;,symsize=0.75  ; 3 = mean ,  4 = min
oplot,listeSA,MsllSA3,psym=8,color=250;,symsize=0.75
oplot,listeSI,MsllSI3,psym=8,color=50;,symsize=0.75
endif 

if (1) then begin
oploterror,listeKOG-0.1,tabstatko(*,0),tabstatko(*,1),psym=8,errcolor=0,errthick=4
oploterror,listeSA+0.1,tabstatsa(*,0),tabstatsa(*,1),psym=8,color=250,errcolor=250,errthick=4
;oploterror,listeSI,tabstatsi(*,0),tabstatsi(*,1),psym=8,color=140,errcolor=140,errthick=4

for i=0,dn-1 do begin

minko=tabstatko(i,2)
minsa=tabstatsa(i,2)

if minko le minsa then begin
oplot,[listeKOG(i)],[minko],psym=2,color=0
listesolkosa(i)=3 ; if ko
endif else begin
oplot,[listeKOG(i)],[minsa],psym=2,color=0
listesolkosa(i)=4 ; if sa
endelse


endfor
save,file="listekosa.sav",listesolkosa,listesolminko,listesolminsa
save,file='whichsol.sav',listesolko,listesolsa,listesolsi
endif

yshift=-10



; Légende
;usersym
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=4;,/fill

oplot,[19.92],[yshift-35],psym=8,symsize=1
xyouts,20.5,yshift-35-0.6,'Init',charsize=1.5

;usersym
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=1,/fill
xyouts,20.5,yshift-39-0.6,'KOGAN',charsize=1.3
oplot,[19.92],[yshift-38.4-0.6],psym=8,symsize=1,color=0
xyouts,20.5,yshift-41-2-0.6,'SA',charsize=1.3
oplot,[19.92],[yshift-42.4-0.6],psym=8,symsize=1,color=250
xyouts,20.5,yshift-41-6-0.6,'Min',charsize=1.3
oplot,[19.92],[yshift-46.4-0.6],psym=2,symsize=1,color=0

device,/close


restore,file="MES-INIT.sav"
tabstatini_FWHM=fltarr(dn,3)
tabstatini_AR=fltarr(dn,3)

for i=0,dn-1 do begin

tabstatini_FWHM(i,0)=mean(TABMESFWHM_INIT(i,*,0))
tabstatini_FWHM(i,1)=sigma(TABMESFWHM_INIT(i,*,0))
tabstatini_FWHM(i,2)=max(TABMESFWHM_INIT(i,*,0))

tabstatini_AR(i,0)=mean(TABMESFWHM_INIT(i,*,5))
tabstatini_AR(i,1)=sigma(TABMESFWHM_INIT(i,*,5))
tabstatini_AR(i,2)=min(TABMESFWHM_INIT(i,*,5))

endfor


set_plot,'ps'
device, yoff=1., ysize=15., xoff=0., xsize=20, /portrait
device,/color,bits_per_pixel=8
device,file='4-ALL-'+'FWHM-'+'5-22.ps'
loadct,39
!p.multi=0
!p.font=0


r=0.7     
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),/fill

plot,listeKOG,tabstatko_FWHM(*,0),title='HPBW = f(Nant) (Kogan/SA)',xtit='Number of antennas',ytit='HPBW ( )',/nodata,yr=[0,90],charsize=2,thick=4,xr=[4,23],/xsty,/ysty  ;!Z(00b0)

oploterror,listeKOG-0.1,tabstatko_FWHM(*,0),tabstatko_FWHM(*,1),psym=8,symsize=1,thick=4
oploterror,listeSA+0.1,tabstatsa_FWHM(*,0),tabstatsa_FWHM(*,1),psym=8,symsize=1,thick=4,color=250,errcolor=250


for i=0,dn-1 do begin
oplot,[listeKOG(i)],[max([tabstatko_FWHM(i,2),tabstatsa_FWHM(i,2)])],psym=2,color=0
endfor

usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=4;,/fill
;oplot,listeKOG,tabstatini_FWHM(*,0),psym=8,symsize=1,thick=4
oploterror,listeKOG,tabstatini_FWHM(*,0),tabstatini_FWHM(*,1),psym=8,symsize=1,thick=4

yshift=100
; Légende
;usersym
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=4;,/fill

oplot,[19.92],[yshift-35],psym=8,symsize=1
xyouts,20.5,yshift-35-0.6,'Init',charsize=1.5

;usersym
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=1,/fill
xyouts,20.5,yshift-39-0.6,'KOGAN',charsize=1.3
oplot,[19.92],[yshift-38.4-0.6],psym=8,symsize=1,color=0
xyouts,20.5,yshift-41-2-0.6,'SA',charsize=1.3
oplot,[19.92],[yshift-42.4-0.6],psym=8,symsize=1,color=250
xyouts,20.5,yshift-41-6-0.6,'Min',charsize=1.3
oplot,[19.92],[yshift-46.4-0.6],psym=2,symsize=1,color=0




usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),/fill
plot,listeKOG,tabstatko_AR(*,0),title='Axial ratio',xtit='Number of antennas',ytit='axial ratio',/nodata,charsize=2,thick=4,xr=[4,23],/xsty,yr=[0.5,2],/ysty
oploterror,listeKOG-0.1,tabstatko_AR(*,0),tabstatko_AR(*,1),psym=8,symsize=1,thick=4
oploterror,listeSA+0.1,tabstatsa_AR(*,0),tabstatsa_AR(*,1),psym=8,symsize=1,thick=4,color=250,errcolor=250

yshift=1.5
; Légende
;usersym
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=4;,/fill

oplot,[19.92],[yshift+0.14],psym=8,symsize=1
xyouts,20.5,yshift+0.14-0.01,'Init',charsize=1.5

;usersym
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=1,/fill
xyouts,20.5,yshift+0.07-0.01,'KOGAN',charsize=1.3
oplot,[19.92],[yshift+0.07],psym=8,symsize=1,color=0
xyouts,20.5,yshift-0.01,'SA',charsize=1.3
oplot,[19.92],[yshift],psym=8,symsize=1,color=250
xyouts,20.5,yshift-0.07-0.01,'Min',charsize=1.3
oplot,[19.92],[yshift-0.07],psym=2,symsize=1,color=0




;for i=0,dn-1 do begin
;oplot,[listeKOG(i)],[min([tabstatko_AR(i,2),tabstatsa_AR(i,2)])],psym=8,color=50
;endfor


usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=4;,/fill
;oplot,listeKOG,tabstatini_AR(*,0),psym=8,symsize=1,thick=4
oploterror,listeKOG,tabstatini_AR(*,0),tabstatini_AR(*,1),psym=8,symsize=1,thick=4

oplot,[0,30],[1,1],line=2,thick=4
device,/close
stop
end
