pro plotbestsol,ENV=ENV
; usage:
Ndeb=5
Nfin=22
Ntry=10
dn=Nfin-Ndeb+1

k=8

; Plot parameters
set_plot,'ps'
!p.font=0
device, yoff=0.8, ysize=32., xoff=0., xsize=22.;, /portrait
device,/color,bits_per_pixel=8
!p.multi=0
!p.font=0

; black point for antenna center 
r=0.5     
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),/fill

dth=2.
tabmesmaxp=fltarr(dn)
tabmesFWHM=fltarr(dn,6)

restore,'whichsol.sav'
restore,'listekosa.sav'

device,file='POWBestsol'+'.ps'
multiplot, [4,5], /square,gap=0.01
for i=0,dn-1 do begin
   
   numstr=strcompress(string(ndeb+i),/remove_all)
   
   tmpa=fltarr(3,ndeb+i)
   

   indkosa=listesolkosa(i)  ; 3 for ko, 4 for sa 

   if indkosa eq 3 then begin
      indtry=listesolminko(i)
      path='progKOG/'
      nameS='saveko-'
      indrec=listesolko(i,indtry)
   endif else if indkosa eq 4 then begin
      indtry=listesolminsa(i)
      path='progSA/'
      nameS='savesa-'
      indrec=listesolsa(i,indtry)
   endif

   trystr=strcompress(string(indtry),/remove_all)
   
   
   restore,file=path+nameS+numstr+'-'+trystr
    
     
   print,'Restoring=',numstr," ",trystr," REC=",indrec eq 1?'END':'RECORD'
   if indrec eq 2 then aa=aarecord    ; 1 for end 2 for rec
      
      ; centering array
      tmpa(0,*)=aa(*,0)-mean(aa(*,0))
      tmpa(1,*)=aa(*,1)-mean(aa(*,1))
      
    maxp=computepowdiag(i,indtry,transpose(reform(tmpa(*,*))),dth,newp,stat,ENV=ENV)
    
    tabmesmaxp(i)=maxp
    tabmesFWHM(i)=stat
     multiplot
  endfor
 
multiplot,/reset
device,/close


device,file='distriBestsol2'+'.ps'
!p.font=0
 multiplot, [4,5], /square,gap=0.01
for i=0,dn-1 do begin
   
   numstr=strcompress(string(ndeb+i),/remove_all)
   
   tmpa=fltarr(3,ndeb+i)
  

   indkosa=listesolkosa(i)  ; 3 for ko, 4 for sa 

   if indkosa eq 3 then begin
      indtry=listesolminko(i)
      path='progKOG/'
      nameS='saveko-'
      indrec=listesolko(i,indtry)
   endif else if indkosa eq 4 then begin
      indtry=listesolminsa(i)
      path='progSA/'
      nameS='savesa-'
      indrec=listesolsa(i,indtry)
   endif

   trystr=strcompress(string(indtry),/remove_all)
   
   restore,file=path+nameS+numstr+'-'+trystr
    
     
   print,'Restoring=',nameS,numstr," ",trystr," REC=",indrec eq 1?'END':'RECORD'
   if indrec eq 2 then aa=aarecord    ; 1 for end 2 for rec
      
         ; centering array
      tmpa(0,*)=aa(*,0)-mean(aa(*,0))
      tmpa(1,*)=aa(*,1)-mean(aa(*,1))
      if i eq 11 then save,file='savesolution16.sav',tmpa
      if i eq 12 then save,file='savesolution17.sav',tmpa
      plot,tmpa(0,*),tmpa(1,*),xr=[-2,2],yr=[-2,2],/xsty,/ysty,psym=8,/isotropic,color=0,background=250,charsize=1.5
      for ii=0,ndeb+i-1 do oplot,tmpa(0,ii)+cos(findgen(181)*2*!dtor)*sqrt(1./!pi/k),tmpa(1,ii)+sin(findgen(181)*2*!dtor)*sqrt(1./!pi/k),color=0
    
      maxp=tabmesmaxp(i)
      xyouts,-1.9,1.55,'SSL='+strcompress(string(maxp,format='(F8.2)'),/remove_all)+' dB',charsize=1.2,charthick=2
      xyouts,-1.9,-1.85,'N='+strcompress(string(ndeb+i,format='(I2)'),/remove_all),charsize=1.5,charthick=2
      
      multiplot
   endfor
   
multiplot,/reset
device,/close




;spawn,'for i in *.ps; do ps2pdf $i;done'

stop
end
