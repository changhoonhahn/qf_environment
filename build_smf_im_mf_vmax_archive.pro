pro build_smf_im_mf_vmax, run, nolit=nolit, twoenv=twoenv
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para        = where(parameters.run eq run)
    param       = parameters[para]
    survey      = strtrim(param.name,1)
    zmin        = param.zmin
    zmax        = param.zmax
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    thresh      = param.thresh
    highenvthresh   = param.highenvthresh
    lowenvthresh    = param.lowenvthresh
    rad_string      = strtrim(string(cylrad),2)
    h_string        = strtrim(string(cylheight),2)
    thresh_string   = strtrim(string(thresh),2)
    
    if (survey EQ 'primus') then begin 
        fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire'] 
        samples = ''
        if keyword_set(nolit) then litsuffix='' $
            else litsuffix = '_lit'
        if keyword_set(nolit) then zbins = mf_zbins(nzbins) $
            else zbins = mf_zbins(nzbins, /literature) 
        if keyword_set(nolit) then zbin = ['0203','0304','0405','05065','06508','0810'] $
            else zbin = ['0204', '0406', '0608', '0810']
        envfname = 'EDP-primus-z0210-numden.fits'
    endif
    if (survey EQ 'sdss') then begin 
        fields = ['sdss']
        samples = '_sdss' 
        zbins = mf_zbins_chh(nzbins, /sdss) 
        zbin = ['nobin']
        litsuffix=''
;       envfname = 'EDP-sdss-z006_0145-numden.fits'
        envfname = 'EDP-sdss-z00375_0145-numden.fits'
    endif 
    struct_print, zbins

    binsize = mf_binsize(bin=0.25,minmass=minmass,maxmass=maxmass)

    sf_q = ['active','quiescent']
    envbin = ['hienv','midenv','lowenv']
    if keyword_set(twoenv) then envbin = ['hienv','lowenv']

    for i=0L, n_elements(sf_q)-1L do begin
        datafile = get_path(/envcount)+'envcount_cylr'+rad_string+'h'+h_string+$
            '_thresh'+thresh_string+samples+'_'+sf_q[i]+litsuffix+'_'+envfname
        data = mrdfits(datafile, 1)
        print, datafile, n_elements(data)
        data = struct_addtags(data,replicate({sfq:''},n_elements(data)))
        data.sfq    = sf_q[i] 
        data.field  = strtrim(data.field,2)
            
        if (i EQ 0L) then all_data = data $
            else all_data = [all_data,data]
    endfor 
    print, 'All Galaxies = ',n_elements(all_data)
    edgecut = where(all_data.edgecut eq 1, nedgecut)
    print, 'All Galaxies within edgecut = ', nedgecut
    all_data = all_data[edgecut]
    
    for iz=0L,nzbins-1L do begin
        zindx   = where((all_data.redshift ge zbins[iz].zlo) and (all_data.redshift lt zbins[iz].zup),zbin_count)
        zbin_data = all_data[zindx]
        print, 'Galaxies within redshift bin '+strmid(strtrim(string(zbins[iz].zlo),2),0,4)+'-'$
            +strmid(strtrim(string(zbins[iz].zup),2),0,4)+' ', zbin_count 
        for ee=0L,n_elements(envbin)-1L do begin
            if keyword_set(twoenv) then begin 
                if (envbin[ee] EQ 'hienv') then env_indx = where(zbin_data.envcount GT highenvthresh, n_envbin)
                if (envbin[ee] EQ 'lowenv') then env_indx = where(zbin_data.envcount LT lowenvthresh, n_envbin) 
            endif else begin
                if ee eq 0 then env_indx = where(zbin_data.envcount GE highenvthresh)
                if ee eq 1 then env_indx = where(zbin_data.envcount GE lowenvthresh AND $ 
                    zbin_data.envcount LT highenvthresh)
                if ee eq 2 then env_indx = where(zbin_data.envcount LT lowenvthresh)
            endelse 
            env_data = zbin_data[env_indx]
            print, 'Galaxies in '+envbin[ee]+' ', n_envbin
            for k=0L,n_elements(sf_q)-1L do begin
                sfqindx = where(env_data.sfq EQ sf_q[k],finaldata_count)
                finaldata = env_data[sfqindx] 
                print, 'Finaldata count = ', finaldata_count
                vmax = finaldata.vmaxavail

                mfdata = im_mf_vmax(finaldata.mass,finaldata.weight/vmax,$
                    binsize=binsize,minmass=minmass,maxmass=maxmass,$
                    masslimit=finaldata.masslimit,rev=rev)
                mfdata = struct_addtags(mfdata,{nfield:intarr(mfdata.nbins),thesefields:strarr(5,mfdata.nbins)})
                print, zbin[iz], sf_q[k], envbin[ee],float(finaldata_count)/float(zbin_count),max(finaldata.envcount),mfdata.ngal

                ; Obtains the fields that contribute to the SMF at each mass bin: 
                if (survey EQ 'sdss') then begin
                    mfdata.nfield = 1
                    mfdata.thesefields[0] = 'sdss'
                endif else begin 
                   for bb=0L,mfdata.nbins-1L do begin 
                       if (rev[bb] ne rev[bb+1]) then begin 
                           allfield = (finaldata.field)[rev[rev[bb]:rev[bb+1]-1]]
                           thesefields = allfield[uniq(allfield,sort(allfield))]
                           mfdata.nfield[bb] = n_elements(thesefields)
                           mfdata.thesefields[0:mfdata.nfield[bb]-1,bb] = thesefields
                        endif  
                    endfor 
                endelse 
                
                ; Since SDSS only has one field, it has to divide SDSS by RA and DEC: 
                if (survey EQ 'sdss') then begin 
                    sdss_hist = hist_nd(transpose([[finaldata.ra],[finaldata.dec]]),$
                        [30D,20D],rev=sdss_rev)
                    notempty = where(sdss_hist gt 100,nfield_jack)
                endif else nfield_jack = n_elements(fields)

                for f=0L,nfield_jack-1L do begin 
                    if (survey EQ 'sdss') then begin 
                        nn = notempty[f] 
                        toss = [sdss_rev[sdss_rev[nn]:sdss_rev[nn+1]-1]]
                        keep = cmset_op(lindgen(n_elements(finaldata)),'and',/not2,toss)
                        reweight = 1D
                    endif else begin 
                        keep = where(strtrim(finaldata.field,2) ne fields[f]) 
                        uniq_field  = (finaldata[keep].field)[uniq(finaldata[keep].field)]
                        uarea       = get_poly_area(/primus,/sr,field=uniq_field)
                        reweight    = get_poly_area(/primus,/sr)/uarea
                    endelse 

                    jack_hist = im_mf_vmax(finaldata[keep].mass,reweight*finaldata[keep].weight/finaldata[keep].vmaxavail,$
                        binsize=binsize,minmass=minmass,maxmass=maxmass,masslimit=finaldata[keep].masslimit)
                    if (f EQ 0L) then begin
                        jack_mfdata = jack_hist
                    endif else begin
                        jack_mfdata = [jack_mfdata, jack_hist]
                    endelse 
                endfor
                
                good = where(mfdata.nfield gt 0, ngood)
                for gg=0L,ngood-1L do begin 
                    if (survey EQ 'sdss') then begin 
                        nbb = nfield_jack
                        m2 = lindgen(nfield_jack)
                    endif else begin 
                        match, strtrim(mfdata.thesefields[*,good[gg]],2),fields ,m1,m2
                        nbb = n_elements(m1)
                    endelse 
                    if (nbb ge 3) then begin 
                        mfdata.phierr_cv[good[gg]] = sqrt((nbb-1.0)/nbb*total((jack_mfdata[m2].phi[good[gg]]$
                            -djs_mean(jack_mfdata[m2].phi[good[gg]]))^2))
                    endif 
                endfor 

                good = where(mfdata.nfield gt 0 and mfdata.phierr_cv gt 0.0) 
                fix = where(mfdata.nfield gt 0 and mfdata.phierr_cv le 0.0,nfix)
                if (nfix ne 0) then begin 
                    get_element, good, fix, nearest
                    mfdata.phierr_cv[fix] = mfdata.phierr_cv[good[nearest]]
                endif 
                
                if (survey EQ 'primus') then fname='smf_primus_cylr'+rad_string+'h'+h_string+$
                    '_thresh'+thresh_string+litsuffix+'_'+strtrim(sf_q[k])+'_'+$
                    envbin[ee]+'_'+zbin[iz]+'z_bin'+string(binsize,format='(F4.2)')+'_'+$
                    strtrim(string(n_elements(envbin)),2)+'envbin.fits'
                if (survey EQ 'sdss') then fname='smf_sdss_cylr'+rad_string+'h'+h_string+$
                    '_thresh'+thresh_string+litsuffix+'_'+strtrim(sf_q[k])+'_'+$
                    envbin[ee]+'_'+zbin[iz]+'z_bin'+string(binsize,format='(F4.2)')+'_'+$
                    strtrim(string(n_elements(envbin)),2)+'envbin.fits'
                print, fname 
                mwrfits, mfdata, get_path(/smf)+fname, /create

                if (sf_q[k] EQ 'active') then begin
                    sf_finaldata = finaldata
                    sf_hist = mfdata 
                endif
                if (sf_q[k] EQ 'quiescent') then begin
                    q_finaldata = finaldata
                    q_hist = mfdata 
                endif
            endfor 

            nbins = sf_hist.nbins
            if (sf_hist.nbins ne q_hist.nbins) then STOP
            qf = {mass:fltarr(nbins),qf:fltarr(nbins),err:fltarr(nbins),err_cv:fltarr(nbins),nfield:intarr(nbins),limit:intarr(nbins),thesefields:strarr(5,nbins)}

            for n=0L,nbins-1L do begin
                qf.mass[n] = q_hist.mass[n]
                qf.qf[n] = q_hist.phi[n]/(q_hist.phi[n]+sf_hist.phi[n])
                qf.limit[n] = min([q_hist.limit[n],sf_hist.limit[n]])
                if (q_hist.phi[n]+sf_hist.phi[n] EQ 0.0) then qf.qf[n] = 0.0

                qf.err[n] = sqrt((sf_hist.phi[n]/(sf_hist.phi[n]+q_hist.phi[n])^2*q_hist.phierr[n])^2$
                    +(q_hist.phi[n]/(sf_hist.phi[n]+q_hist.phi[n])^2*sf_hist.phierr[n])^2)
                if (q_hist.phi[n]+sf_hist.phi[n] EQ 0.0) then qf.err[n] = 0.0

                nfieldcount = 0L
                for nn=0L,4L do begin 
                    if (q_hist.thesefields[nn,n] NE sf_hist.thesefields[nn,n]) then begin
                        qf.thesefields[nn,n] = strtrim(q_hist.thesefields[nn,n]+sf_hist.thesefields[nn,n],2)
                    endif else begin 
                        qf.thesefields[nn,n] = strtrim(q_hist.thesefields[nn,n],2)
                    endelse 
                    if qf.thesefields[nn,n] ne '' then nfieldcount = nfieldcount+1
                endfor 
                qf.nfield[n] = nfieldcount
            endfor 
                
            if (survey EQ 'sdss') then begin 
                sf_sdss_hist = hist_nd(transpose([[sf_finaldata.ra],[sf_finaldata.dec]]),$
                    [30D,20D],rev=sf_rev)
                q_sdss_hist = hist_nd(transpose([[q_finaldata.ra],[q_finaldata.dec]]),$
                    [30D,20D],rev=q_rev)
                sf_nbins= n_elements(sf_sdss_hist)
                q_nbins = n_elements(q_sdss_hist)
                notempty = where(sf_sdss_hist gt 100 and q_sdss_hist gt 100,nfield_jack)
            endif else nfield_jack = n_elements(fields)

            for f=0L,nfield_jack-1L do begin 
                if (survey EQ 'sdss') then begin 
                    nn = notempty[f] 
                    sf_toss = [sf_rev[sf_rev[nn]:sf_rev[nn+1]-1]]
                    q_toss = [q_rev[q_rev[nn]:q_rev[nn+1]-1]]
                    sf_keep = cmset_op(lindgen(n_elements(sf_finaldata)),'and',/not2,sf_toss)
                    q_keep = cmset_op(lindgen(n_elements(q_finaldata)),'and',/not2,q_toss)
                    sf_reweight = 1D
                    q_reweight = 1D
                endif else begin 
                    sf_keep = where(strtrim(sf_finaldata.field,2) ne fields[f]) 
                    q_keep  = where(strtrim(q_finaldata.field,2) ne fields[f]) 
                    sf_uniq_field = (sf_finaldata[sf_keep].field)[uniq(sf_finaldata[sf_keep].field)]
                    q_uniq_field  = (q_finaldata[q_keep].field)[uniq(q_finaldata[q_keep].field)]
                    sf_uarea = get_poly_area(/primus,/sr,field=sf_uniq_field)
                    q_uarea  = get_poly_area(/primus,/sr,field=q_uniq_field)
                    sf_reweight = get_poly_area(/primus,/sr)/sf_uarea
                    q_reweight  = get_poly_area(/primus,/sr)/q_uarea
                endelse 
                sf_jack_hist = im_mf_vmax(sf_finaldata[sf_keep].mass,sf_reweight*sf_finaldata[sf_keep].weight/sf_finaldata[sf_keep].vmaxavail,$
                    binsize=binsize,minmass=minmass,maxmass=maxmass,masslimit=sf_finaldata[sf_keep].masslimit)
                q_jack_hist = im_mf_vmax(q_finaldata[q_keep].mass,q_reweight*q_finaldata[q_keep].weight/q_finaldata[q_keep].vmaxavail,$
                    binsize=binsize,minmass=minmass,maxmass=maxmass,masslimit=q_finaldata[q_keep].masslimit)
                jack_qf = {qf:fltarr(nbins)}
                for n=0L,nbins-1L do begin 
                    jack_qf.qf[n] = q_jack_hist.phi[n]/(q_jack_hist.phi[n]+sf_jack_hist.phi[n])
                    if (q_jack_hist.phi[n]+sf_jack_hist.phi[n] EQ 0.0) then jack_qf.qf[n] = 0.0
                endfor 

                if (f EQ 0L) then begin
                    jack_qfdata = jack_qf
                endif else begin
                    jack_qfdata = [jack_qfdata, jack_qf]
                endelse 
            endfor

            good = where(qf.nfield gt 0, ngood) 
            for gg=0L,ngood-1L do begin 
                if (survey EQ 'sdss') then begin 
                    nbb = nfield_jack
                    m2 = lindgen(nfield_jack)
                endif else begin 
                    match, strtrim(qf.thesefields[*,good[gg]],2), fields, m1, m2
                    nbb = n_elements(m2)
                endelse 
                if (nbb GE 3) then begin 
                    qf.err_cv[good[gg]] = sqrt((nbb-1.0)/nbb*$
                        total((jack_qfdata[m2].qf[good[gg]]-$
                        djs_mean(jack_qfdata[m2].qf[good[gg]]))^2))
                endif 
            endfor 
            good = where(qf.nfield gt 0 and qf.err_cv gt 0.0) 
            fix = where(qf.nfield gt 0 and qf.err_cv le 0.0, nfix)
            if (nfix ne 0) then begin 
                get_element, good, fix, nearest
                qf.err_cv[fix] = qf.err_cv[good[nearest]]
            endif 
           
            if (survey EQ 'primus') then qf_fname='qf_primus_cylr'+rad_string+'h'+h_string+$
                '_thresh'+thresh_string+litsuffix+'_'+envbin[ee]+'_'+zbin[iz]+$
                'z_bin'+string(binsize,format='(F4.2)')+'_'+strtrim(string(n_elements(envbin)),2)+'envbin.fits'
            if (survey EQ 'sdss') then qf_fname='qf_sdss_cylr'+rad_string+'h'+h_string+$
                '_thresh'+thresh_string+litsuffix+'_'+envbin[ee]+'_'+zbin[iz]+$
                'z_bin'+string(binsize,format='(F4.2)')+'_'+strtrim(string(n_elements(envbin)),2)+'envbin.fits'
            print, qf_fname
            mwrfits, qf, get_path(/smf)+qf_fname, /create
        endfor
    endfor 
end

;    if keyword_set(prank) then begin 
;        lowthresh = 20
;        print, '     z_low      z_high      zero percent'
;        for k=0L,nzbins-1L do begin
;            edgecut = allgal_combine.edgecut eq 1
;            zindx   = allgal_combine.redshift ge zbins[k].zlo and allgal_combine.redshift lt zbins[k].zup
;            zdata   = allgal_combine[where(edgecut and zindx, zcount)]
;
;            zeroenv = zdata.envcount eq 0 
;            zeroenvdata = zdata[where(zeroenv, zerocount)]
;
;            zeropercent = 100.0*float(zerocount)/float(zcount)
;            print, zbins[k].zlo, zbins[k].zup, zeropercent
;            lowthresh = lowthresh>zeropercent
;        endfor
;        print, lowthresh
;        lowthresh = float(ceil(lowthresh))
;        print, 'Low environment percentage rank threshold=', lowthresh
;    endif 

;        if keyword_set(prank) then begin
;            zindx_data = allgal_combine[where(edgecut and zindx, zindx_count)]
;            zindx_data.prank = get_percentage_rank(zindx_data.envcount)
;
;            for j=0L,n_elements(sf_q)-1L do begin
;                sfqindx = zindx_data.sfq eq sf_q[j]
;                
;                ;if keyword_set(literature) then minmass = get_min_masslimit(iii,sf_q[j],/literature) $
;                ;    else minmass = get_min_masslimit(iii,sf_q[j])
;                smcomp = zindx_data.mass gt zindx_data.masslimit
;
;                for jj=0L,n_elements(hilow)-1L do begin 
;                    lowthresh = 20.0>min(zindx_data.prank)
;;                    print, lowthresh
;                    if jj eq 0 then prankindx = zindx_data.prank gt 80.0
;                    if jj eq 1 then prankindx = zindx_data.prank gt lowthresh and zindx_data.prank le 80.0
;                    if jj eq 2 then prankindx = zindx_data.prank le lowthresh
;
;                    prank_data = zindx_data[where(sfqindx and smcomp and prankindx, prankcount)]
;                    print, sf_q[j],zbin[iii],hilow[jj],100.0*float(prankcount)/float(zindx_count),' percent, min=',min(prank_data.envcount),',max=',max(prank_data.envcount)
;                    
;                    hist = im_mf_vmax(prank_data.mass,prank_data.weight/prank_data.vmaxavail,masslimit=prank_data.masslimit)
;                    
;                    if (survey EQ 'sdss') then fname='smf_primus_cylr'+radius_string+'h'+height_string+$
;                        '_thresh'+threshold_string+'_nbin'+nbin_string+litsuffix+strtrim(sf_q[j])+'_prank_'+$
;                        hilow[jj]+'_'+zbin[iii]+'z.fits'
;                    if (survey EQ 'sdss') then fname='smf_sdss_cylr'+radius_string+'h'+height_string+$
;                        '_thresh'+threshold_string+'_nbin'+nbin_string+litsuffix+strtrim(sf_q[j])+'_prank_'+$
;                        hilow[jj]+'_'+zbin[iii]+'z.fits'
;                    mwrfits, hist, get_path(/smf)+fname, /create
;                endfor
;            endfor
;        endif else begin       
