pro prof_plot
;+
; Name
;      prof_plot
; Purpose
;      Plot the velocity and temperature profiles of the
;      adiabatic expansion motion, and compare it to the
;      exponential approximation.
;
; Written by sw, Dec. 04, 2018
;-

  ;System parameter
  rootpath='~/data/emcal/model/cc85/analytical/adiabatic/'+ $
    'adiabatic_nei/adiabatic/model_check/'
  figpath =rootpath+'figures/'
  
  ;Plot parameter
  set_plot,'ps'
  device,/encapsulated,preview=2,bits=8,/color
  tvlct,[0,255,  0,  0,255,255,  0], $
        [0,  0,255,  0,  0,127,255], $
        [0,  0,  0,255,255,  0,255]
  ;   black, R,  G,  B,meg,org,cyan
  charsize=1.5 & notesize=1.2 & charthick=3
  xthick=3 & ythick=3 & thick=3
  xtitle='Radius (pc)' & ytitle='Ratio'
  
  ;Physical parameters defined in Ji+2006, Sec.3.4
  mu=1.4
  u0=1e3 ;km/s
  T0=5e6 ;K
  cs0=sqrt(5*!k*T0/(3*mu*!mp))/1e5 ;km/s
  Tc=T0+mu*!mp*(u0*1e5)^2/(5*!k) ;K
  r0=0.3 ;pc
  mdot0=3e-5*!msolar/!yr ;g/s
  c1=mdot0/(4*!dpi) ;constant coefficient for \rho*r^2*u
  nH0=c1/(u0*1e5*(r0*!pc)^2*mu*!mp) ;in cm-3
  
  ;Derive the radial profile
  u =sqrt(u0^2+3*cs0^2-3*cs0^2/(findgen(200)/5+1)^(4/3.)) ;in km/s
  ur=r0/((1+(u0^2-u^2)/3/cs0^2)^3*(u/u0)^2)^0.25 ;in pc
  ; print,u[0:103]
  ; print,1/(u[0:103]/u0)^2
  ; print,1/(1+(u0^2-u[0:103]^2)/3/cs0^2)^3
  rho=c1/(u*1e5*(ur*!pc)^2) ;in g/cm3
  nH=rho/(mu*!mp) ;in cm-3
  ; T =10^(-findgen(221)/100.)*T0 ;in K
  T =(reverse(dindgen(400)+1)/400.)^(4/3.)*T0 ;in K
  Tr=r0/((T/T0)^3*(Tc-T)/(Tc-T0))^0.25 ;in pc
  nT=n_elements(T)
  
  ;Deal with the dynamical timescale
  nu=n_elements(u)
  ur_u=[ur[1:(nu-1)],ur[nu-1]]
  del_ur=ur_u-ur
  del_t =del_ur/u*!pc/1e5
  del_tau=del_t*nH
  
  xrange=[0.95,300]*r0 & yrange=[0.85,1.1] & resrange=[-1,1]*0.12
  u_app =sqrt(3*cs0^2+u0^2-3*cs0^2/(ur/r0)^(4/3.))
  nH_app=u0/u_app*(r0/ur)^2*nH0
  T_app =T0*(r0/Tr)^(4/3.)
  device,filename=figpath+'adia.exp_phy.eps'
  plot,xrange,yrange,/nodata,xrange=xrange,yrange=yrange, $
    xstyle=1,ystyle=1,title='!6',xtitle='',ytitle=ytitle, $
    ymargin=[0,0.5],xthick=xthick,ythick=ythick, $
    xtickname=replicate(' ',2),charsize=charsize, $
    charthick=charthick,/xlog,ylog=0, $
    position=[0.15,0.3,0.95,0.95]
  oplot,ur,u/u0,thick=thick,color=1
  oplot,ur,u_app/u0,thick=thick*2,linestyle=2,color=1
  oplot,ur,nH/nH0*(ur/r0)^2,thick=thick,color=2
  oplot,ur,nH_app/nH0*(ur/r0)^2,thick=thick*2,linestyle=2,color=2
  oplot,Tr,T/T0*(Tr/r0)^(4/3.),thick=thick,color=3
  oplot,Tr,T_app/T0*(Tr/r0)^(4/3.),thick=thick*2,linestyle=2,color=3
  oplot,[0,0]+r0,yrange,linestyle=2,thick=thick
  legend,['!8v!6/!8v!6!I0!N', $
    '(!8n!I!6H!N/!8n!6!IH0!N)(!8r!6/!8r!6!I0!N)!U2!N', $
    '(!8T!6/!8T!6!I0!N)(!8r!6/!8r!6!I0!N)!U4/3!N' ], $
    thick=[3,3,3],color=[1,2,3],linestyle=[0,0,0], $
    charsize=notesize,charthick=charthick,bthick=xthick,/bottom
  plot,xrange,[0,0],linestyle=2,thick=thick,xrange=xrange, $
    yrange=resrange,xstyle=1,ystyle=1,xtitle=xtitle,ytitle='res', $
    xthick=xthick,ythick=ythick,yminor=2,yticks=2,ytickv=[-0.1,0,0.1], $
    xticklen=!P.TICKLEN*4.3,charsize=charsize, $
    charthick=charthick,/noerase,position=[0.15,0.15,0.95,0.30],/xlog
  oplot,[0,0]+r0,resrange,linestyle=2,thick=thick
  oplot,ur,(u-u_app)/u,thick=thick,color=1
  oplot,ur,(nH-nH_app)/nH,thick=thick,color=2
  oplot,Tr,(T-T0*(r0/Tr)^(4/3.))/T,thick=thick,color=3
  device,/close
  
  ;Log the physical parameters
  u_app =sqrt(3*cs0^2+u0^2-3*cs0^2/(Tr/r0)^(4/3.))
  nH_app=u0/u_app*(r0/Tr)^2*nH0
  openw,lun,rootpath+'adia.exp_phy.info',/get_lun
  printf,lun,'R          velo       dens       kT'
  for i=0L,nT-1 do printf,lun,Tr[i]*!pc,u_app[i]*1e5, $
    nH_app[i],T[i]*!k/!keV,format='(e10.4,3(1x,e10.4))'
  free_lun,lun
  
  device,filename=figpath+'delta_r.eps'
  Tr_u=[Tr[1:(nT-1)],Tr[nT-1]]
  plot,Tr,Tr_u-Tr,xrange=xrange,yrange=[1e-3,1e1],xstyle=1,ystyle=1, $
    /xlog,/ylog,xtitle=xtitle,ytitle='!7D!8r!6 (pc)',charsize=charsize, $
    charthick=charthick,xthick=xthick,ythick=ythick,thick=thick
  ; oplot,Tr,(Tr_u[0]-Tr[0])/Tr[0]*Tr,thick=3,color=1,linestyle=2
  ; oplot,Tr,(Tr_u[1]-Tr[1])*(Tr[1]/Tr)^2,thick=3,color=1,linestyle=2
  ; legend,['Real Case','!7D!8r!M?!8r!6!U-2!N'],color=[0,1],linestyle=[0,2], $
  ;   thick=[0,0]+thick,charsize=notesize,charthick=charthick;,/right
  device,/close
  
  tau_app=(Tr_u-Tr)*!pc/u_app/1e5*nH_app
  device,filename=figpath+'tau_val.eps'
  plot,Tr,tau_app/1e7,xrange=xrange,yrange=[1.5,3.5],xstyle=1,ystyle=1, $
    xtitle=xtitle,ytitle='!8n!6!IH!N!MX!7D!8r!6/!8v!6 (10!U7!N cm!U-3!Ns)', $
    /xlog,charsize=charsize,charthick=charthick,xthick=xthick, $
    ythick=ythick,thick=thick
  ; oplot,Tr,(Tr_u[0]-Tr[0])/Tr[0]*Tr,thick=3,color=1,linestyle=2
  ; oplot,Tr,(Tr_u[1]-Tr[1])*(Tr[1]/Tr)^2,thick=3,color=1,linestyle=2
  ; legend,['Real Case','!7D!8r!M?!8r!6!U-2!N'],color=[0,1],linestyle=[0,2], $
  ;   thick=[0,0]+thick,charsize=notesize,charthick=charthick;,/right
  device,/close
  
  accum_tau_app=total(tau_app,/cumulative)
  device,filename=figpath+'accum_tau_val.eps'
  plot,Tr,accum_tau_app/1e8,xrange=xrange,yrange=[1.,100.],xstyle=1,ystyle=1, $
    title='!6',xtitle=xtitle,/xlog,/ylog,charsize=charsize,charthick=charthick, $
    xthick=xthick,ythick=ythick,thick=thick, $
    ytitle='Cumulative !8n!6!IH!N!MX!7D!8r!6/!8v!6 (10!U8!N cm!U-3!N s)'
  device,/close
  
  
  set_plot,'x'

END