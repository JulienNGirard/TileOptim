pro plotres,REC=REC
; usage:
Ndeb=5
Nfin=22
Ntry=10
dn=Nfin-Ndeb+1
liste=indgen(dn)+ndeb

namerec=""
if REC then namerec="REC-"



; Loading savesets
restore,file='MES-'+namerec+'KOG.sav'
restore,file='BestKOGSLL.sav'
MsllKO=tabmesmaxp
restore,file='MES-'+namerec+'SA.sav'
restore,file='BestSASLL.sav'
MsllSA=tabmesmaxp
restore,file='MES-'+namerec+'SI.sav'
restore,file='BestSISLL.sav'
MsllSI=tabmesmaxp

; loading masks
if REC then begin
maskKOG=tabKOGmaskrec
maskSA=tabSAmaskrec
maskSI=tabSImaskrec
endif else begin
maskKOG=tabKOGmaskend
maskSA=tabSAmaskend
maskSI=tabSImaskend
endelse

MsslKO2=fltarr(dn)
MsslSA2=fltarr(dn)
MsslSI2=fltarr(dn)

for i=0,dn-1 do begin
MsslKO2(i)=min(MsllKO(i,maskKOG))
MsslSA2(i)=min(MsllSA(i,maskSA))
MsslSI2(i)=min(MsllSI(i,maskSI))
endfor


set_plot,'ps'
device, yoff=1., ysize=15., xoff=0., xsize=20, /portrait
device,/color,bits_per_pixel=8
device,file='3-ALL-'+namerec+'5-22.ps'
loadct,39
!p.multi=0
!p.font=1





; usersym
r=1
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=4 ;,/fill

; Plot initial distribution SLL
plot,listeKOG,tabKOGSLL(0,*),title='Lowest SLL = f(Nant)',xtit='Number of antennas',ytit='SLL (dB)',/nodata,yr=[-70,max(tabSASLL(0,*))],charsize=2,thick=4,xr=[4,23],/xsty
oploterror,listeKOG,tabKOGSLL(0,*),tabKOGSLL(1,*),psym=8,symsize=1,thick=4

; Fit a curve on it
yfit=linfit(listeKOG,reform(tabKOGSLL(0,*)))
print,'Fit coefficient',' slope=',yfit(1),' y0=',yfit(0)
oplot,listeKOG,yfit(1)*listeKOG+yfit(0),thick=4,line=1

; usersym
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),thick=1,/fill


oplot,listeKOG,MsslKO2,psym=8;,symsize=0.75  ; 3 = mean ,  4 = min
oplot,listeSA,MsslSA2,psym=8,color=250;,symsize=0.75
;oplot,listeSI,MsslSI2,psym=8,color=50;,symsize=0.75

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
;xyouts,20.5,yshift-41-6-0.6,'Simple',charsize=1.3
;oplot,[19.92],[yshift-46.4-0.6],psym=8,symsize=1,color=50

device,/close



stop
end
