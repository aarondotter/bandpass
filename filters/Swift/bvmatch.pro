PRO kuruczuvot,modfile,avin,vmatch,bvmatch

readcol,modfile,garbage,p,lambda,nu,flux,format='a,i,f,f,f,',skipline=2

lambda=lambda*10

flux=4*flux*(3e18)/lambda^2

readcol,'uvw2.may09.dat',w2l,w2f
readcol,'uvm2.may09.dat',m2l,m2f
readcol,'uvw1.oct09.dat',w1l,w1f
uur = MRDFITS('swuuu_20041120v104.arf',1)
bbr = MRDFITS('swubb_20041120v104.arf',1)
vvr = MRDFITS('swuvv_20041120v104.arf',1)

uul = uur.WAVE_MAX + (uur.WAVE_MIN-uur.WAVE_MAX)/2.0
bbl = bbr.WAVE_MAX + (bbr.WAVE_MIN-bbr.WAVE_MAX)/2.0
vvl = vvr.WAVE_MAX + (vvr.WAVE_MIN-vvr.WAVE_MAX)/2.0

;Create dust models in step of 0.05 mag up to 1.50
av = FINDGEN(31)*0.05
nsp=n_elements(flux)
d=fltarr(31,nsp)
duse=fltarr(nsp)

done=0

while done eq 0 do begin
   peidust,1,avin,lambda/10000,duse

;dust up the model
;f=flux*d[5,*]
f=flux*duse
;now interpolate the model to the points of the filter function
spw2 = INTERPOL(f,lambda,w2l)
spm2 = INTERPOL(f,lambda,m2l)
spw1 = INTERPOL(f,lambda,w1l)
spuu = INTERPOL(f,lambda,uul)
spbb = INTERPOL(f,lambda,bbl)
spvv = INTERPOL(f,lambda,vvl)

;plot,uul,uur.SPECRESP,xrange=[1700,6000]
;oplot,bbl,bbr.SPECRESP
;oplot,vvl,vvr.SPECRESP
;oplot,w2l,w2f
;oplot,m2l,m2f
;oplot,w1l,w1f

abw2=-2.5*ALOG10(TOTAL(10.*w2l*spw2*w2f)/$
                 TOTAL(10.*w2f/w2l))
abm2=-2.5*ALOG10(TOTAL(10.*m2l*spm2*m2f)/$
                 TOTAL(10.*m2f/m2l))
abw1=-2.5*ALOG10(TOTAL(10.*w1l*spw1*w1f)/$
                 TOTAL(10.*w1f/w1l))
abuu=-2.5*ALOG10(TOTAL((uur.WAVE_MIN-uur.WAVE_MAX)*uul*spuu*uur.SPECRESP)/$
                 TOTAL((uur.WAVE_MIN-uur.WAVE_MAX)*uur.SPECRESP/uul))
abbb=-2.5*ALOG10(TOTAL((bbr.WAVE_MIN-bbr.WAVE_MAX)*bbl*spbb*bbr.SPECRESP)/$
                 TOTAL((bbr.WAVE_MIN-bbr.WAVE_MAX)*bbr.SPECRESP/bbl))
abvv=-2.5*ALOG10(TOTAL((vvr.WAVE_MIN-vvr.WAVE_MAX)*vvl*spvv*vvr.SPECRESP)/$
                 TOTAL((vvr.WAVE_MIN-vvr.WAVE_MAX)*vvr.SPECRESP/vvl))

; mab=mvega+mab of vega
; so vega mags = ab mags - vega
vvv=abvv+0.01
vbb=abbb+0.12
vuu=abuu-1.02
vw1=abw1-1.48
vm2=abm2-1.71
vw2=abw2-1.72

distance=vmatch-vvv

print,'W2=',vw2+distance
print,'M2=',vm2+distance
print,'W1=',vw1+distance
print,'U=',vuu+distance
print,'B=',vbb+distance
print,'V=',vvv+distance
print,''
print,'B-V=',vbb-vvv
print,'M2-U=',vm2-vuu
print,avin/3.1

coldiff=bvmatch-vbb+vvv
if abs(coldiff) gt 0.01 then avin=avin+coldiff*3.1
if abs(coldiff) lt 0.01 then done=1
print,coldiff

change=''
read,change
if change eq '1' then done=1

endwhile

end
