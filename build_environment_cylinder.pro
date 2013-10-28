pro build_environment_cylinder,run,Nransack,primus=primus, sdss=sdss, literature=literature
    if keyword_set(primus) then sample = ['']
    if keyword_set(primus) then fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
    if keyword_set(sdss) then sample = ['_sdss']
    if keyword_set(sdss) then fields = ['sdss']
    subsample = ['all', 'active', 'quiescent']

    t0 = systime(1)
    print, t0
    
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para = where(parameters.run eq run)
    param = parameters[para]
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    nbin        = param.nbin
    threshold   = param.thresh
    thrshld     = float(param.thresh)/100.0

    if keyword_set(primus) then zbins = mf_zbins_chh(nzbins, zmin=zmin, zmax=zmax, /literature)
    if keyword_set(sdss) then zbins = mf_zbins_chh(nzbins, zmin=zmin, zmax=zmax, /sdss)
    
    radius_string = strtrim(string(cylrad),2)
    height_string = strtrim(string(cylheight),2)
    threshold_string = strtrim(string(threshold),2)
    nbin_string = strtrim(string(nbin),2)

    if keyword_set(primus) then envfname='EDP-primus-z0210-numden.fits' 
    if keyword_set(sdss) then envfname='EDP-sdss-z006_0145-numden.fits' 
    envfile = get_path(/envt)+envfname
    print, envfile
    print, 'zmin=',zmin,' zmax=',zmax,' cylrad=',cylrad,' cylheight=',cylheight,' nbin=',nbin,' threshold=',thrshld

    survey=''
    for i=0L,n_elements(fields)-1L do survey = survey+fields[i]+'_'
    ransackfile     = 'ransack_'+survey+strtrim(string(Nransack),2)+'.fits'
    print, ransackfile
    ransack_data    = mrdfits(get_path(/ransack)+ransackfile,1) 
    ran_ra          = ransack_data.ra
    ran_dec         = ransack_data.dec

    for i=0L,n_elements(sample)-1L do begin
        for ii = 0L,n_elements(subsample)-1L do begin
            if keyword_set(literature) then ext='_lit' else ext=''
            mffile = get_path(/target)+'target_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
                '_nbin'+nbin_string+'_'+subsample[ii]+sample[i]+ext+'.fits'
            print, mffile
            mfdata = mrdfits(mffile, 1)
            ngal = n_elements(mfdata)
            target = replicate({ra:0.D, dec:0.D, zprimus:0.}, ngal)

            target.ra       = mfdata.ra
            target.dec      = mfdata.dec
            target.zprimus  = mfdata.redshift 

            nmatch = replicate({envcount:0.}, n_elements(target))
            nmatch.envcount = get_environment_cylinder(envfile=envfile, targ=target, rad=cylrad, h=cylheight)
            
            if keyword_set(primus) then edgecut = get_edgecut(ran_ra, ran_dec, target ,zmin=zmin,zmax=zmax,$
                rad=cylrad,nbin=nbin,thresh=thrshld,/primus)
            if keyword_set(sdss) then edgecut = get_edgecut(ran_ra, ran_dec, target ,zmin=zmin,zmax=zmax,$
                rad=cylrad,nbin=nbin,thresh=thrshld,/sdss)
            
            print, 'time=', (systime(1)-t0)/60.0
       
            output  = struct_addtags(struct_addtags(mfdata, nmatch), edgecut)
            fname   = get_path(/envcount)+'envcount_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+'_nbin'$
                +nbin_string+sample[i]+'_'+subsample[ii]+ext+'_'+envfname
            print, fname
            mwrfits, output, fname, /create
        endfor
   endfor
   print, 'total time=', (systime(1)-t0)/60.0
end 
