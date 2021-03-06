pro build_envcount_qf, run, nolit=nolit, twoenv=twoenv
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
    endif
    if (survey EQ 'sdss') then begin 
        fields = ['sdss']
        samples = '_sdss' 
        zbins = mf_zbins_chh(nzbins, /sdss)
        zbin = ['nobin']
        litsuffix=''
    endif 

    binsize = mf_binsize(bin=0.25,minmass=minmass,maxmass=maxmass)

    sf_q = ['active','quiescent']
    if keyword_set(twoenv) then envbin = ['hienv','lowenv'] $
        else envbin = ['hienv','midenv','lowenv']
    
    for iz=0L,nzbins-1L do begin
        for ee=0L,n_elements(envbin)-1L do begin
; import environment SMF data from build_environment_smf.pro
            smfdata_sf_file = get_path(/envsmfdata)+'smf_'+survey+'_cylr'+rad_string+'h'+h_string+$
                '_thresh'+thresh_string+litsuffix+'_active_'+$
                envbin[ee]+'_'+zbin[iz]+'z_'+strtrim(string(n_elements(envbin)),2)+'envbin.fits' 
            smfdata_q_file = get_path(/envsmfdata)+'smf_'+survey+'_cylr'+rad_string+'h'+h_string+$
                '_thresh'+thresh_string+litsuffix+'_quiescent_'+$
                envbin[ee]+'_'+zbin[iz]+'z_'+strtrim(string(n_elements(envbin)),2)+'envbin.fits' 
            sf_finaldata = mrdfits(smfdata_sf_file,1)
            q_finaldata = mrdfits(smfdata_q_file,1)
            print, smfdata_sf_file
            print, smfdata_q_file

; import SMF files generated by build_envcount_smf.pro
            smf_sf_file='smf_'+survey+'_cylr'+rad_string+'h'+h_string+$
                '_thresh'+thresh_string+litsuffix+'_active_'+$
                envbin[ee]+'_'+zbin[iz]+'z_bin'+string(binsize,format='(F4.2)')+'_'+$
                strtrim(string(n_elements(envbin)),2)+'envbin.fits'
            smf_q_file='smf_'+survey+'_cylr'+rad_string+'h'+h_string+$
                '_thresh'+thresh_string+litsuffix+'_quiescent_'+$
                envbin[ee]+'_'+zbin[iz]+'z_bin'+string(binsize,format='(F4.2)')+'_'+$
                strtrim(string(n_elements(envbin)),2)+'envbin.fits'
            sf_hist = mrdfits(get_path(/smf)+smf_sf_file,1) 
            q_hist = mrdfits(get_path(/smf)+smf_q_file,1) 
            print, smf_sf_file 
            print, smf_q_file 

            nbins = sf_hist.nbins
            if (sf_hist.nbins ne q_hist.nbins) then STOP ; SF and Q smf files don't match up!
            if (array_equal(sf_hist.mass, q_hist.mass) NE 1) then STOP ; double checking!

            qf = {mass:fltarr(nbins),qf:fltarr(nbins),err:fltarr(nbins),err_cv:fltarr(nbins),nfield:intarr(nbins),limit:intarr(nbins),thesefields:strarr(5,nbins)}

            for n=0L,nbins-1L do begin
                qf.mass[n] = q_hist.mass[n]
                qf.qf[n] = q_hist.phi[n]/(q_hist.phi[n]+sf_hist.phi[n])
                if (q_hist.phi[n]+sf_hist.phi[n] EQ 0.0) then qf.qf[n] = 0.0

                qf.limit[n] = max([q_hist.limit[n],sf_hist.limit[n]])   ; mass limit is the max mass limit

                qf.err[n] = sqrt((sf_hist.phi[n]/(sf_hist.phi[n]+q_hist.phi[n])^2*q_hist.phierr[n])^2$
                    +(q_hist.phi[n]/(sf_hist.phi[n]+q_hist.phi[n])^2*sf_hist.phierr[n])^2)
                if (q_hist.phi[n]+sf_hist.phi[n] EQ 0.0) then qf.err[n] = 0.0

                q_thesefields = strtrim(q_hist.thesefields[*,n],2)
                q_n0fields = q_thesefields[where(q_thesefields NE '', n_qnotempty)]
                sf_thesefields = strtrim(sf_hist.thesefields[*,n],2)
                sf_n0fields = sf_thesefields[where(sf_thesefields NE '', n_sfnotempty)]

                thesefields = [q_n0fields, sf_n0fields]
                ;print, '----------------------------------'
                ;print, thesefields
                uniq_thesefields = thesefields[uniq(thesefields[sort(thesefields)])] 
                if (n_elements(uniq_thesefields) EQ 5L) then qf.thesefields[*,n] = uniq_thesefields $
                    else qf.thesefields[*,n] = [ uniq_thesefields, strarr(5L-n_elements(uniq_thesefields)) ] 
                if ((n_qnotempty EQ 0) OR (n_sfnotempty EQ 0)) then qf.nfield[n] = 0 $
                    else qf.nfield[n] = n_elements(uniq_thesefields)
                
                ;print, 'q', q_hist.thesefields[*,n]
                ;print, 'sf', sf_hist.thesefields[*,n]
                ;print, 'q and sf', cmset_op(strtrim(q_hist.thesefields[*,n],2),'and',strtrim(sf_hist.thesefields[*,n],2))
                ;print, 'qf', qf.thesefields[*,n]
            endfor 
                
            if (survey EQ 'sdss') then begin 
                min_ra = min([min(sf_finaldata.ra), min(q_finaldata.ra)])
                min_dec = min([min(sf_finaldata.dec), min(q_finaldata.dec)])
                max_ra = max([max(sf_finaldata.ra), max(q_finaldata.ra)])
                max_dec = max([max(sf_finaldata.dec), max(q_finaldata.dec)])

                sf_sdss_hist = hist_nd(transpose([[sf_finaldata.ra],[sf_finaldata.dec]]),$
                    [30D,20D],min=[min_ra,min_dec],max=[max_ra,max_dec],rev=sf_rev)
                q_sdss_hist = hist_nd(transpose([[q_finaldata.ra],[q_finaldata.dec]]),$
                    [30D,20D],min=[min_ra,min_dec],max=[max_ra,max_dec],rev=q_rev)
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
                    sf_uarea = get_poly_area(/primus,/sr,field=strtrim(sf_uniq_field,2))
                    q_uarea  = get_poly_area(/primus,/sr,field=strtrim(q_uniq_field,2))
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
            fixx = where(qf.nfield gt 0 and qf.err_cv le 0.0, nfix)
            if (nfix ne 0) then begin 
                get_element, good, fixx, nearest
                qf.err_cv[fixx] = qf.err_cv[good[nearest]]
            endif 
           
            qf_fname='qf_'+survey+'_cylr'+rad_string+'h'+h_string+$
                '_thresh'+thresh_string+litsuffix+'_'+envbin[ee]+'_'+zbin[iz]+$
                'z_bin'+string(binsize,format='(F4.2)')+'_'+strtrim(string(n_elements(envbin)),2)+'envbin.fits'
            print, qf_fname
            mwrfits, qf, get_path(/qfenvcount)+qf_fname, /create
        endfor
    endfor 
end
