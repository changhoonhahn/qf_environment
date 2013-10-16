;+
; NAME:
;   get_edgecut
;
; PURPOSE:
;   Determines whether or not the input points are on the edges of the polygons. First, ransack is used 
;   to generate random points on the 2mask polygons. Then spherematch is used to determine how many of the 
;   random points are within a specified distance (mlength(z), a function of redshift) of the input points. 
;   The input objects are classified into bins according to redshift, with each bin having an average mlength(z).
;
; CALLING SEQUENCE:
;   get_edgecount_bin_SMF(input=, zmin=, zmax=, nbin=, thresh=)
;
; INPUTS:
;   input - data structure input
;   zmin - Minimum redshift value of interest.
;   zmax - Maximum redshift value of interest.
;   nbin - Number of redshift bins we want to divide the data.
;
; OUTPUTS:
;   nmatch - number of environment objects inside the cylinder centered around each of the target objects
;
; FUNCTIONS USED:
;   rd_tfile('data.fits',3)  [Requires str2arr.pro and deriv_arr.pro]
;   propmotdis(z, OmegaM, OmegaL, weq=weq)
;   get_total_poly_area() 
;
; PROCEDURES USED:
;   spherematch, ra1, dec1, ra2, dec2, matchlength, match1, match2, distance12, maxmatch=
;
;-------------------------------------------------------

function get_edgecut,ran_ra,ran_dec, target,zmin=zmin,zmax=zmax,rad=rad,nbin=nbin,thresh=thresh,primus=primus,sdss=sdss
    if (n_elements(ran_ra) eq 0 OR n_elements(ran_dec) EQ 0 OR n_elements(target) EQ 0) then print, 'INPUT ERROR'

    ransack_num = n_elements(ran_ra)
    print, 'Nransack=',ransack_num

;Arranging the data into redshift bins:
    zdiff=zmax-zmin
    zstep=zdiff/nbin
    b=dblarr(nbin+1)
    b[nbin]=zmax
    z0=zmin

    for j=0L,nbin-1 do begin
    b[j]=z0
    z0=z0+zstep
    endfor

    outputs=replicate({edgecut:0L}, n_elements(target))
    binstart=0
    bincount=0

; Length at different redshift for each of the bins
; Note, the average of the redshift bins are used as redshift: 
    mlength=dblarr(nbin)
    for k=0L, nbin-1 do begin
        mlength[k]=(float(rad)/(3000.0*propmotdis(0.5*(b[k]+b[k+1]), 0.3, 0.7)))*(180.0/!PI)
    endfor

;Total area of all the polygons:
    if keyword_set(primus) then totalarea=get_poly_area(/primus)
    if keyword_set(sdss) then totalarea=get_poly_area(/sdss)

    for i=0L,nbin-1 do begin
;        help, /mem
    	zbin = where(target.zprimus ge b[i] and target.zprimus lt b[i+1], zbincount)
     	bin = target[zbin]
        print, 'bin dimensions=',zbincount 

        binra   = bin.ra
        bindec  = bin.dec
        binz    = bin.zprimus
     	spherematch, ran_ra, ran_dec, binra, bindec, mlength[i], m_ran, m_targ,$
            bindist12, maxmatch=0
     
;     	ml      = (float(rad)/(3000.0*propmotdis(binz[m_targ], 0.3, 0.7)))*(180.0/!PI)
;     	ikeep   = where(bindist12 lt ml, nkeep)
        print, 'bindist=', size(bindist12, /dimensions);, 'nkeep=', nkeep

;        mt      = m_targ[ikeep]
        mt      = m_targ
        isort   = sort(mt)
        sorted  = mt[sort(mt)]
        iuniq   = uniq(mt[isort])
        istart  = 0L
        nmatch  = lonarr(zbincount)
 
        for m=0L, n_elements(iuniq)-1L do begin
            iend    = iuniq[m]
            icurr   = isort[istart:iend]
            nmatch[mt[icurr[0]]] = n_elements(icurr)
            istart  = iend+1L
        endfor 

        threshold = (float(ransack_num)/totalarea)*(mlength[i]^2*!PI)*thresh
; We determine the (number of random points within mlength)/(expected number of random points): 
;        randense = (float(nmatch)/(mlength[i]^2*!PI))/(float(ransack_num)/totalarea)
; We determine a cut for whether the point is on the edge or not on the edge: 
        edgecut=lonarr(zbincount)
        edgecut[where(nmatch lt threshold)] = 1
;        for n=0L, zbincount-1L do begin
;            if nmatch[n] lt threshold[i] then begin
;                    edgecut[n]=0
;            endif else begin
;                    edgecut[n]=1
;            endelse
;        endfor

        print, 'average nmatch=', mean(nmatch),' threshold=',threshold,' mean edgecut=',mean(edgecut)

        binend      = binstart+zbincount-1
        outputs[binstart:binend].edgecut = edgecut
        binstart    = binstart+zbincount
        bincount    = bincount+zbincount
    endfor
    return, outputs
end 
