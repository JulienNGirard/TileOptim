pro reduce,ENV=ENV,MODE=MODE
;usage: reduce,/env,MODE='KOG'
Ndeb=5
Nfin=22
Ntry=10
dn=Nfin-Ndeb+1
liste=indgen(dn)+ndeb

tabSLLg=fltarr(8,dn)
tabitnum=fltarr(4,dn)
tabmaskrec=fltarr(Ntry,dn)
tabmaskend=fltarr(Ntry,dn)

case MODE of   
   'KOG':  begin
      nameP='savePKOG-'
      nameS='saveko-'
   end
   'SA':    begin
      nameP='savePSA-'
      nameS='savesa-'
   end
   'SI':   begin
      nameP='savePSI-'
      nameS='savesi-'
   end
   else:    begin
      print,'Error in MODE - stopping'
      stop
   end

endcase


for ni=0,dn-1 do begin
   
   restore,file=nameP+strcompress(string(Ndeb+ni),/remove_all)
   
   
      tabpmaxrec=tabtabout(7,*)    ; Valeur record
      tabpmaxend=tabtabout(1,*)    ; Valeur de fin

; filtrage valeurs abérrantes
   if ENV then begin
      maskend=where(tabpmaxend lt -10 and tabpmaxend gt -150)
      maskrec=where(tabpmaxrec lt -10 and tabpmaxrec gt -150)
 
      if total(maskend) eq -1 then maskend=where(tabpmaxend lt -10)
      if total(maskrec) eq -1 then maskrec=where(tabpmaxrec lt -10)
   endif else begin
      maskend=where(tabpmaxend lt 0)
      maskrec=where(tabpmaxrec lt 0)
   endelse

   tabmaskend(maskend,ni)=1
   tabmaskrec(maskrec,ni)=1

   tabSLLg(0,ni)=mean(tabtabout(2,*))  ; valeur de départ moyenne
   tabSLLg(1,ni)=sigma(tabtabout(2,*)) ; sigma
   
   tabSLLg(2,ni)=mean(tabpmaxend(maskend))  ; valeur SLL END moyenne
   tabSLLg(3,ni)=sigma(tabpmaxend(maskend)) ; sigma
   tabSLLg(4,ni)=min(tabpmaxend(maskend))   ; min
   
   tabSLLg(5,ni)=mean(tabpmaxrec(maskrec))  ; valeur SLL REC moyenne
   tabSLLg(6,ni)=sigma(tabpmaxrec(maskrec)) ; sigma
   tabSLLg(7,ni)=min(tabpmaxrec(maskrec))   ; min

   
   for ntry=0,9 do begin                    ; sur chaque essai
      restore,file=nameS+strcompress(string(Ndeb+ni),/remove_all)+'-'+strcompress(string(ntry),/remove_all)
      maxprecord=min(tabsll)                ; record SLL
      tmploc=where(tabsll eq maxprecord)    ; localisation record
      
      if ENV then begin
         if maxprecord lt -10 then begin
            if ntry eq 0 then tablocrecord=tmploc else tablocrecord=[tablocrecord,tmploc]
         endif
      endif else begin
         if maxprecord lt 0 then begin
            if ntry eq 0 then tablocrecord=tmploc else tablocrecord=[tablocrecord,tmploc]
         endif
      endelse

   endfor

   tabitnum(0,ni)=mean(tablocrecord)    ; valeur moyenne du record
   tabitnum(1,ni)=sigma(tablocrecord)   ; sigma
   tabitnum(2,ni)=mean(tabtabout(8,maskend)) ; valeur moyenne du nombre final d'itérations
   tabitnum(3,ni)=sigma(tabtabout(8,maskend)) ; sigma
endfor

; Sauvegarde

case MODE of   
   'KOG':  begin
      
      listeKOG=liste            ; liste des N_ant
      tabKOGSLL=tabSLLg         ; valeurs SLL
      tabKOGitnum=tabitnum
      tabKOGmaskend=tabmaskend
      tabKOGmaskrec=tabmaskrec
      save,file='Best'+MODE+'SLL.sav',tabKOGSLL,listeKOG,tabKOGmaskend,tabKOGmaskrec,tabKOGitnum

   end
   'SA':    begin
      listeSA=liste            ; liste des N_ant
      tabSASLL=tabSLLg          ; valeurs SLL
      tabSAitnum=tabitnum
      tabSAmaskend=tabmaskend
      tabSAmaskrec=tabmaskrec
      save,file='Best'+MODE+'SLL.sav',tabSASLL,listeSA,tabSAmaskend,tabSAmaskrec,tabSAitnum

   end
   'SI':   begin
      listeSI=liste            ; liste des N_ant
      tabSISLL=tabSLLg          ; valeurs SLL
      tabSIitnum=tabitnum
      tabSImaskend=tabmaskend
      tabSImaskrec=tabmaskrec
      save,file='Best'+MODE+'SLL.sav',tabSISLL,listeSI,tabSImaskend,tabSImaskrec,tabSIitnum

   end
   else:    begin
      print,'Error in MODE - stopping'
      stop
   end   
endcase



set_plot,'ps'
device, yoff=1., ysize=15., xoff=0., xsize=20, /portrait
device,/color,bits_per_pixel=8
device,file='1-'+MODE+'-5-22-SLL.ps'
loadct,39
!p.multi=0
!p.font=0

plot,liste,tabSLLg(0,*),title=MODE+' - SLL = f(Nant)',xtit='Number of antennas',ytit='SLL (dBc)',/nodata,yr=[min(tabSLLg(2,*)),max(tabSLLg(0,*))],charsize=2,thick=4,xr=[4,23],/xsty

oploterror,liste,tabSLLg(0,*),tabSLLg(1,*),psym=6,symsize=0.75

yfit=linfit(liste,reform(tabSLLg(0,*)))
print,yfit
oplot,liste,yfit(1)*liste+yfit(0),thick=4,line=2
oploterror,liste,tabSLLg(2,*),tabSLLg(3,*),psym=5,symsize=0.75

oplot,[20],[-3],psym=6
oplot,[20],[-5],psym=5
xyouts,20.5,-3.5,'Start SLL'
xyouts,20.5,-5.6,'Best SLL'
device,/close

stop
end
