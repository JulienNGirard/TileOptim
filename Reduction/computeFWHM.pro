function computeFWHM,pdb,dth

  x1=indgen(90./dth+1)
  x2=findgen(10000)*(90./dth+1)/10000
            

  stat=fltarr(6)
 ; tabmaxFWHM=fltarr(nd,6)
 ; tabaxialratio=fltarr(nd,2)
            tabFWHMslice=fltarr(360./dth+1)
            for iphi=0,360./dth do begin
               
               pdbslice=pdb(*,iphi)
               pdbint=interpol(pdbslice,x1,x2)
               
               mask=where(pdbint gt -3.)
               
               if total(mask) ne -1 then begin
          
                  nm=n_elements(mask)
                  for imask=0,nm-2 do begin
                  diff=mask(imask+1)-mask(imask)
                  cnt=imask
                  if abs(diff) gt 1 then break
                  endfor

                  fwhm=dth*x2(cnt) ; in degrees   
            endif else begin
               fwhm=-42
            endelse
            tabFWHMslice(iphi)=2*fwhm

            endfor
            w=where(tabFWHMslice ge 0)
            stat(0)=mean(tabFWHMslice(w))
            stat(1)=sigma(tabFWHMslice(w))
            stat(2)=min(tabFWHMslice(w))
            stat(3)=max(tabFWHMslice(w))
            stat(4)=stat(3)-stat(2)
            stat(5)=stat(3)/stat(2)
         
         
         ;locbest=(where(tabFWHM(i,masktry,4) eq min(tabFWHM(i,masktry,4))))(0)
         ;locbest2=(where(tabFWHM(i,masktry,3) eq max(tabFWHM(i,masktry,3))))(0)
         
         ;tabmaxFWHM(i,0)=tabFWHM(i,locbest,0)
         ;tabmaxFWHM(i,1)=tabFWHM(i,locbest,1)
         ;tabmaxFWHM(i,2)=tabFWHM(i,locbest,2)
         ;tabmaxFWHM(i,3)=tabFWHM(i,locbest,3)
         ;tabmaxFWHM(i,4)=tabFWHM(i,locbest,4)

        ; tabaxialratio(i,0)=mean(tabFWHM(i,masktry,5))
         ;tabaxialratio(i,1)=sigma(tabFWHM(i,masktry,5))
         
;save,file='FWHM'+MODE,tabFWHM,tabmaxFWHM,tabaxialratio

return,stat

end
