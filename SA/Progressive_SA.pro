Pro Progressive_SA,ndep,n,P
r=2
k=8
distmin=2*sqrt(1./!pi/k)

restore,'Init-distrib.sav'
tabpmax=fltarr(P)
tabimin=fltarr(P)

  for i=Ndep,Ndep+n-1 do begin
     tabadeb=fltarr(3,i,P)
     tabaend=fltarr(3,i,P)
     tabtabout=fltarr(9,P)
     nom="Confsa_"+strcompress(string(long(i)),/remove_all)
     na=i

     for numtry=0,P-1 do begin
        a=tabini(*,0:i-1,i-5,numtry)
        tabadeb(*,*,numtry)=a
        array_beam_sa_repuls2,numtry,na,a,k,aout,tabout,Xcos=2;,/plot
        print,'Nant=',i,'Ntry=',numtry
        tabPmax(numtry)=tabout(1)
        tabimin(numtry)=tabout(0)
        tabAend(*,*,numtry)=aout
        tabtabout(*,numtry)=tabout
     endfor
        save,file='savePSA-'+strcompress(string(na),/remove_all),na,tabpmax,tabimin,tabadeb,tabaend,tabtabout
        print,'test LOOP '+strcompress(string(i)+'/'+string(100.*(i-3)/(n-3)),/remove_all)+' % DONE'
  endfor
print,"Study successfully done"
  stop
end
