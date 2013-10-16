function mf_select_sfr_quiescent_chh, sfr, mass, z=z, active=active
; jm12mar21ucsd - select quiescent and actively star-forming galaxies
;   using an evolving cut in SFR/M

    mfpath = '/mount/moon1/ioannis/research/projects/primus/mf/2165/'
    bimod = mrdfits(mfpath+'sfr_bimodality.fits.gz',1,/silent)

    sfrperpcut = poly(z-bimod.zpivot,bimod.linecoeff)

    sfrperp = sfr-bimod.slope*(mass-bimod.masspivot)
    iquiescent = where(sfrperp lt sfrperpcut,comp=iactive)
    if keyword_set(active) then indx = iactive else $
      indx = iquiescent
return, indx
end
