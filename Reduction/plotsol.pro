pro plotsol,REC=REC,ENV=ENV,MODE=MODE
; usage:
Ndeb=5
Nfin=22
Ntry=10
dn=Nfin-Ndeb+1

k=8

; Plot parameters
set_plot,'ps'
device, yoff=0.8, ysize=32., xoff=0., xsize=22.;, /portrait
device,/color,bits_per_pixel=8
!p.multi=0
!p.font=1

; black point for antenna center 
r=0.5     
usersym,r*sin(findgen(17)*!pi/8),r*cos(findgen(17)*!pi/8),/fill

if REC then begin
namefile="distriREC-"
namefilepow="Pow-REC-"
namerec="REC-"
endif else begin
namefile="distri-"
namefilepow="Pow-"
namerec=""
endelse

case MODE of   
   'KOG':  begin
      nameP='savePKOG-'
      nameS='saveko-'
     ; restore,file="Best"+MODE+"SLL.sav"
     ; tabSLLg=tabKOGSLL
   end
   'SA':    begin
      nameP='savePSA-'
      nameS='savesa-'
      ;restore,file="Best"+MODE+"SLL.sav"
      ;tabSLLg=tabSASLL
   end
   'SI':   begin
      nameP='savePSI-'
      nameS='savesi-'
     ; restore,file="Best"+MODE+"SLL.sav"
     ; tabSLLg=tabSISLL
   end
   else:    begin
      print,'Error in MODE - stopping'
      stop
   end

endcase


dth=2.
tabmesmaxp=fltarr(dn,ntry)
tabmesFWHM=fltarr(dn,ntry,6)
tabmesAR=fltarr(dn,2)
for i=0,dn-1 do begin
   
   numstr=strcompress(string(ndeb+i),/remove_all)
   device,file=namefilepow+numstr+'.ps'
   
   tmpa=fltarr(3,ndeb+i,ntry)
   multiplot, [3,4], /square,gap=0.01
   restore,file=nameP+numstr
   
   for itry=0,ntry-1 do begin
      trystr=strcompress(string(itry),/remove_all)
      print,numstr," ",trystr
      restore,file=nameS+numstr+'-'+trystr
      if REC then aa=aarecord      

      
      ; centering array
      tmpa(0,*,itry)=aa(*,0)-mean(aa(*,0))
      tmpa(1,*,itry)=aa(*,1)-mean(aa(*,1))
      
    maxp=computepowdiag(i,itry,transpose(reform(tmpa(*,*,itry))),dth,newp,stat,ENV=ENV)
    
    tabmesmaxp(i,itry)=maxp
    tabmesFWHM(i,itry,*)=stat
     multiplot
  endfor

   
multiplot,/reset
endfor

device,/close
save,file='MES-'+namerec+MODE+'.sav',tabmesmaxp,tabmesFWHM


for i=0,dn-1 do begin
   
   numstr=strcompress(string(ndeb+i),/remove_all)
   device,file=namefile+numstr+'.ps'
   
   tmpa=fltarr(3,ndeb+i,ntry)
   multiplot, [3,4], /square    ;,gap=0.01
   restore,file=nameP+numstr
   
   for itry=0,ntry-1 do begin
      trystr=strcompress(string(itry),/remove_all)
      print,numstr," ",trystr
      restore,file=nameS+numstr+'-'+trystr
      if REC then aa=aarecord
         ;maxp=tabtabout(7,itry)  ; min maxp record
     ; endif else begin
     ;    maxp=tabtabout(1,itry)  ; min maxp end
     ; endelse
      
      ; centering array
      tmpa(0,*,itry)=aa(*,0)-mean(aa(*,0))
      tmpa(1,*,itry)=aa(*,1)-mean(aa(*,1))
      
      plot,tmpa(0,*,itry),tmpa(1,*,itry),xr=[-2,2],yr=[-2,2],/xsty,/ysty,psym=8,/isotropic,color=0,background=250,charsize=1.5
      for ii=0,ndeb+i-1 do oplot,tmpa(0,ii,itry)+cos(findgen(181)*2*!dtor)*sqrt(1./!pi/k),tmpa(1,ii,itry)+sin(findgen(181)*2*!dtor)*sqrt(1./!pi/k),color=0
    
      maxp=tabmesmaxp(i,itry)
      xyouts,-1.9,1.65,'SSL='+strcompress(string(maxp,format='(F8.2)'),/remove_all)+' dB',charsize=1,charthick=2
      
      multiplot
   endfor
   
multiplot,/reset
endfor





spawn,'for i in *.ps; do ps2pdf $i;done'

stop
end
