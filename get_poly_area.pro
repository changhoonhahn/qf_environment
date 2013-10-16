;+
; NAME:
;   get_poly_area
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; OUTPUTS:
;
; PROCEDURES USED:
;   read_mangle_polygons, infile, polygons, id [, unit=]
;
; MODIFICATION HISTORY:
;--------------------------------------------------------
function get_poly_area, primus=primus, sdss=sdss, sr=sr, field=field
    if keyword_set(primus) then begin 
        dir ='/global/data/scr/chh327/primus/survey_regions/fields/'
        fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
        if keyword_set(field) then fields = [field]
        areas=dblarr(n_elements(fields))

        for i=0L,n_elements(fields)-1L do begin
            read_mangle_polygons, dir+fields[i]+'_field_galex_window_2mask.ply',polygon
            areas[i]=total(polygon.str,/double)
        endfor

        totarea=total(areas)
        if keyword_set(sr) then begin
            return, totarea
        endif else begin
            return, (totarea)*(180.0/!PI)^2
        endelse
    endif

    if keyword_set(sdss) then begin
;        dir = '/mount/moon1/ioannis/research/projects/primus/mf/2165/'
        dir = '/global/data/scr/chh327/primus/science/mf/2165/'
        read_mangle_polygons, dir+'dr72bsafe0_galex_final.ply', sdss
        area = total(sdss.str, /double)
        if keyword_set(sr) then begin
            return, area
        endif else begin
            return, area*(180.0/!PI)^2
        endelse 
    endif
end

