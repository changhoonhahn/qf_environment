function get_vmax_avail, z, vz
    vmaxavail = interpol(vz.v_max_avail,vz.z,z,/spline)
    return, vmaxavail
end 
pro build_target_sample,run,Nransack,Nrandom,literature=literature, primus=primus, sdss=sdss
    mfdatapath  = get_path(/mfdata)
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para = where(parameters.run eq run)
    param= parameters[para]
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    threshold   = param.thresh
    nbin        = param.nbin
    h100 = mf_h100()
    if keyword_set(primus) then begin 
        zbins = mf_zbins_chh(nzbins, zmin=zmin, zmax=zmax, /literature) 
        jm_zbins = mf_zbins(jm_nzbins,zmin=jm_zmin,zmax=jm_zmax, /literature)
    endif
    if keyword_set(sdss) then begin
        zbins = mf_zbins_chh(nzbins, zmin=zmin, zmax=zmax, /sdss) 
        jm_zbins = mf_zbins(jm_nzbins,zmin=jm_zmin,zmax=jm_zmax, /sdss)
    endif 

    struct_print, zbins
    struct_print, jm_zbins

    radius_string   = strtrim(string(cylrad),2)
    height_string   = strtrim(string(cylheight),2)
    threshold_string= strtrim(string(threshold),2)
    nbin_string     = strtrim(string(nbin),2)
    rsk_string      = strtrim(string(Nransack),1)
    ran_string      = strtrim(string(Nrandom),1)
    
    if keyword_set(primus) then area=get_poly_area(/primus,/sr)
    if keyword_set(sdss) then area=get_poly_area(/sdss,/sr)

    if keyword_set(sdss) then sample = ['_sdss']
    if keyword_set(sdss) then fields = ['sdss']
    if keyword_set(primus) then sample = ['']
    if keyword_set(primus) then fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire'] 
    subsample   = ['all', 'active', 'quiescent']

    survey=''
    for i=0L,n_elements(fields)-1L do begin
        survey=survey+fields[i]+'_'
    endfor

;Vmax,avail file generated from build_vmax_avail.pro
    vmaxavail_name=get_path(/vmaxavail)+'vmax_avail_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+'_nbin'+$
        nbin_string+'_'+survey+'ran'+ran_string+'_rsk'+rsk_string+'.fits'
    print, vmaxavail_name
    vz = mrdfits(vmaxavail_name,1)

    for i=0L,n_elements(sample)-1L do begin
        for ii =0L, n_elements(subsample)-1L do begin 
            if keyword_set(literature) then ext = '_lit' else ext = '' 
            mfdatafile = mfdatapath+'mfdata_'+subsample[ii]+'_supergrid01'+sample[i]+ext+'.fits.gz'
            print, mfdatafile
            mfdata = mrdfits(mfdatafile,1)

            zlim = mfdata.z GE zmin AND mfdata.z LT zmax
            data = mfdata[where(zlim)]
            ngal = n_elements(data)	
            data = struct_addtags(data,replicate({vmax:0.,vmaxavail:0.},ngal))

            for j=0L,nzbins-1L do begin 
                zindx = where(data.z GT zbins[j].zlo AND data.z LE zbins[j].zup,zindx_count)
                data_zbin = data[zindx]

                jm_zmin_evol = data_zbin.zmin_evol>jm_zbins[j].zlo
                jm_zmax_evol = data_zbin.zmax_evol<jm_zbins[j].zup
                zmin_evol = data[zindx].zmin_evol>zbins[j].zlo
                zmax_evol = data[zindx].zmax_evol<zbins[j].zup
                vmax_evol = data[zindx].vmax_evol
;                data[zindx].vmax        = (vmax_evol/(lf_comvol(jm_zmax_evol)-lf_comvol(jm_zmin_evol)))$
;                    *(lf_comvol(zmax_evol)-lf_comvol(zmin_evol))
                data[zindx].vmax = area/3.0*(lf_comvol(zmax_evol)-lf_comvol(zmin_evol))*(1/h100)^3.0
                data[zindx].vmaxavail   = (vmax_evol/(lf_comvol(jm_zmax_evol)-lf_comvol(jm_zmin_evol))) $
                    *(get_vmax_avail(zmax_evol,vz)-get_vmax_avail(zmin_evol,vz))
                if (min(data[zindx].vmaxavail) LT 0.0) then begin
                    print, "ERROR Negative Vmax,avail"
                    print, "Try using fewer points for Vmax,avail interpol"
                    STOP
                endif 
                lowmasstest_indx = where(data.z GT zbins[j].zlo AND data.z LE zbins[j].zup AND data.mass LT 9.5, lowmasstest_count)
                lowmasstest = data[lowmasstest_indx]
;                for jj=0L,lowmasstest_count-1L do begin 
;                    print, lowmasstest[jj].mass
;                    print, jm_zbins[j].zlo, jm_zbins[j].zup
;                    print, zbins[j].zlo, zbins[j].zup
;                    print, lowmasstest[jj].zmin_evol, lowmasstest[jj].zmax_evol
;                    print, lowmasstest[jj].vmax_evol, lowmasstest[jj].vmax
;                endfor 
            endfor 
           
            target = replicate({ra:0.D,dec:0.D,redshift:0.,mr_01:0.,mg_01:0.,mb_00:0.,mass:0.,masslimit:0.,age:0.,SFR:0.,$ 
                       field:' ',weight:0.,vmax_evol:0.,vmax:0.,vmaxavail:0.},ngal)
            if keyword_set(sdss) then target.field = 'sdss' else target.field = data.field 
            target.ra       = data.ra
            target.dec      = data.dec
            target.redshift = data.z
            target.mr_01    = data.mr_01
            target.mg_01    = data.mr_01+data.gmr_01
            target.mb_00    = data.mb_00
            target.mass     = data.mass
            target.masslimit= data.masslimit
            target.age      = data.age
            target.SFR      = data.SFR
            target.weight   = data.weight
            target.vmax_evol= data.vmax_evol
            target.vmax     = data.vmax
            target.vmaxavail= data.vmaxavail
            
            dir = get_path(/target)
            if keyword_set(literature) then begin
                fname = dir+'target_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
                    '_nbin'+nbin_string+'_'+subsample[ii]+sample[i]+'_lit.fits'
            endif else begin
                fname = dir+'target_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
                    '_nbin'+nbin_string+'_'+subsample[ii]+sample[i]+'.fits'
            endelse
            print, fname
            mwrfits, target, fname, /create
        endfor 
    endfor 
end
