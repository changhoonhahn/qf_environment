;This code will generate the table for V_max,avail. 
pro build_vmax_avail,run,n_rsk,n_ran,primus=primus,sdss=sdss
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para = where(parameters.run eq run)
    param = parameters[para]
    zmin = param.zmin 
    zmax = param.zmax
    nbin = param.nbin
    cylrad = param.cylrad
    cylheight = param.cylheight
    threshold = param.thresh
    thrshld = float(param.thresh)/100.0
    nz = 50 
    
    radius_string       = strtrim(string(cylrad),2)
    height_string       = strtrim(string(cylheight),2)
    threshold_string    = strtrim(string(threshold),2)
    nbin_string         = strtrim(string(nbin),2)
    ran_string          = strtrim(string(n_ran),1)
    rsk_string          = strtrim(string(n_rsk),1)

    rsk_dir = get_path(/ransack)
    if keyword_set(primus) then fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
    if keyword_set(sdss) then fields = ['sdss']

    survey = ''
    for i=0L,n_elements(fields)-1L do begin 
        survey = survey+fields[i]+'_'
    endfor
    ransack_fname   = rsk_dir+'ransack_'+survey+strtrim(string(n_rsk),1)+'.fits'
    random_fname    = rsk_dir+'random_'+survey+strtrim(string(n_ran),1)+'.fits'
  
    print, ransack_fname
    print, random_fname
    ransack = mrdfits(ransack_fname,1)
    random  = mrdfits(random_fname,1)

    ran_vol = randomn(systime(1),n_ran,/uniform) 
    ran_z   = dblarr(n_ran)

    vmax_vmin   = lf_comvol(zmax)-lf_comvol(zmin)
    vmin        = lf_comvol(zmin)
    for i=0L,n_ran-1L do begin  
        ran_vol[i]  = ran_vol[i]*(vmax_vmin)+vmin
        ran_z[i]    = lf_vtoz(ran_vol[i])
    endfor    
    randomdata = replicate({ra:0.,dec:0.,zprimus:0.},n_ran) 
    randomdata.ra       = random.ra
    randomdata.dec      = random.dec
    randomdata.zprimus  = ran_z

    rsk_ra = ransack.ra
    rsk_dec= ransack.dec

    if keyword_set(primus) then edgecut = get_edgecut(rsk_ra,rsk_dec,randomdata,zmin=zmin,zmax=zmax,$
        rad=cylrad,nbin=nbin,thresh=thrshld,/primus)
    if keyword_set(sdss) then edgecut = get_edgecut(rsk_ra,rsk_dec,randomdata,zmin=zmin,zmax=zmax,$
        rad=cylrad,nbin=nbin,thresh=thrshld,/sdss)
    ran_data = struct_addtags(randomdata,edgecut)

    zbin     = dblarr(nz)
    eff_vmax = dblarr(nz)
    for i=0L,nz-1L do begin 
        zbin[i] = zmin+float(i+1L)*(zmax-zmin)/float(nz)
        z_indx  = ran_data.zprimus LE zbin[i] AND ran_data.zprimus GE zmin
        edgecut_indx = ran_data.edgecut eq 1L

        edgecut_pts     = ran_data[where(z_indx AND edgecut_indx, edgecut_count)] 
        unedgecut_pts   = ran_data[where(z_indx, unedgecut_count)]
        eff_vmax[i] = float(edgecut_count)/float(unedgecut_count)*lf_comvol(zbin[i])
        print, zbin[i],eff_vmax[i]/lf_comvol(zbin[i]),float(edgecut_count)/float(unedgecut_count)
    endfor 
    out = replicate({z:0.,v_max_avail:0.},nz)
    out.z = zbin
    out.v_max_avail = eff_vmax

    outname = get_path(/vmaxavail)+'vmax_avail_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
        '_nbin'+nbin_string+'_'+survey+'ran'+ran_string+'_rsk'+rsk_string+'.fits'
    print, outname
    mwrfits,out,outname,/create
end 
