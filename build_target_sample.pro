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
    if keyword_set(primus) then zbins = mf_zbins_chh(nzbins, zmin=zmin, zmax=zmax, /literature) 
    if keyword_set(sdss) then zbins = mf_zbins_chh(nzbins, zmin=zmin, zmax=zmax, /sdss) 

    radius_string   = strtrim(string(cylrad),2)
    height_string   = strtrim(string(cylheight),2)
    threshold_string= strtrim(string(threshold),2)
    nbin_string     = strtrim(string(nbin),2)
    rsk_string      = strtrim(string(Nransack),1)
    ran_string      = strtrim(string(Nrandom),1)
    
    if keyword_set(primus) then area=get_poly_area(/primus,/sr)
    if keyword_set(sdss) then area=get_poly_area(/sdss,/sr)

    if keyword_set(sdss) then sample = ['sdss']
    if keyword_set(primus) then sample = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire'] 
    subsample   = ['all', 'active', 'quiescent']

    survey=''
    for i=0L,n_elements(sample)-1L do begin
        survey=survey+sample[i]+'_'
    endfor

;Vmax,avail file generated from build_vmax_avail.pro
    vmaxavail_name=get_path(/vmaxavail)+'vmax_avail_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+'_nbin'+$
        nbin_string+'_'+survey+'ran'+ran_string+'_rsk'+rsk_string+'.fits'
    print, vmaxavail_name
    vz = mrdfits(vmaxavail_name,1)

    for i=0L,n_elements(sample)-1L do begin
        for ii =0L, n_elements(subsample)-1L do begin 
            if keyword_set(literature) then ext = '_lit' else ext = '' 
            mfdatafile = mfdatapath+'mfdata_'+subsample[ii]+'_supergrid01_'+sample[i]+ext+'.fits.gz'
            mfdata = mrdfits(mfdatafile,1)

            zlim = mfdata.z GE zmin AND mfdata.z LT zmax
            data = mfdata[where(zlim)]
            ngal = n_elements(data)	
           
            target = replicate({class:' ',ra:0.D,dec:0.D,redshift:0.,mr_01:0.,mb_00:0.,mass:0.,masslimit:0.,age:0.,SFR:0.,$ 
                       field:' ',weight:0.,vmax:0.,vmaxavail:0.},ngal)
            target.field    = sample[i]
            target.ra       = data.ra
            target.dec      = data.dec
            target.redshift = data.z
            target.mr_01    = data.mr_01
            target.mb_00    = data.mb_00
            target.mass     = data.mass
            target.masslimit= data.masslimit
            target.age      = data.age
            target.SFR      = data.SFR
            target.weight   = data.weight
            target.vmax     = data.vmax_evol
           
            vmax_avail = dblarr(ngal)
            for j=0L,nzbins-1L do begin 
                zindx = data.z GT zbins[j].zlo AND data.z LE zbins[j].zup
                data_zbin = data[where(zindx,zindx_count)]

                zlo = zbins[j].zlo
                zup = zbins[j].zup
                zmin_evol = data_zbin.zmin_evol
                zmax_evol = data_zbin.zmax_evol
                vmax_evol = data_zbin.vmax_evol
                vmax_avail = (vmax_evol/(lf_comvol(zmax_evol)-lf_comvol(zmin_evol)))$
                *(get_vmax_avail(zmax_evol<zup,vz)-get_vmax_avail(zmin_evol>zlo,vz))

                target[where(target.redshift GT zbins[j].zlo AND target.redshift LE zbins[j].zup,target_zbin_count)].vmaxavail = vmax_avail
            endfor 

;            for j=0L,ngal-1 do begin
;                if keyword_set(sdss) then print, j
;                if keyword_set(sdss) then zbin = get_zbin((data.z)[j],/sdss)
;                if keyword_set(primus) then zbin = get_zbin((data.z)[j],/primus)
;                zlo = zbin[0]
;                zhi = zbin[1] 
;                zmin_evol = (data.zmin_evol)[j]
;                zmax_evol = (data.zmax_evol)[j]
;                vmax_evol = (data.vmax_evol)[j]
;                vmax[j] = (vmax_evol/(lf_comvol(zmax_evol)-lf_comvol(zmin_evol)))$
;                *(get_vmax_avail(zmax_evol<zhi,vz)-get_vmax_avail(zmin_evol>zlo,vz))
;                if (vmax[j] LT 0.0) then print, "ERROR"
;                if keyword_set(sdss) then vmax[j]=(area/3.0)*(get_vmax_avail(zmax_evol<zhi,/sdss)$
;                    -get_vmax_avail(zmin_evol>zlo,/sdss))*(1.0/h100)^3.0
;                if keyword_set(primus) then vmax[j]=(area/3.0)*(get_vmax_avail(zmax_evol<zhi,/primus)$
;                    -get_vmax_avail(zmin_evol>zlo,/primus))*(1.0/h100)^3.0
;                    if (zmin_evol lt zlo) or (zmax_evol gt zhi) then begin
;                        print, 'Vmax changed', target[j].redshift, zmin_evol, zlo, zmax_evol, zhi
;                        vmax[j] = (vmax_evol/(lf_comvol(zmax_evol)-lf_comvol(zmin_evol)))$
;                                *(lf_comvol(zmax_evol<zhi)-lf_comvol(zmin_evol>zlo))
;                    endif else begin
;                        vmax[j] = vmax_evol
;                    endelse
;            endfor

            dir = get_path(/target)
            if keyword_set(literature) then begin
                fname = dir+'target_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
                    '_nbin'+nbin_string+'_'+subsample[ii]+'_'+sample[i]+'_lit_test.fits'
            endif else begin
                fname = dir+'target_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
                    '_nbin'+nbin_string+'_'+subsample[ii]+'_'+sample[i]+'_test.fits'
            endelse
            mwrfits, target, fname, /create
        endfor 
    endfor 
end
