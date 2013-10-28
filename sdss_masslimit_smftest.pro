pro sdss_masslimit_smftest
    mfdata_dir = "/global/data/scr/chh327/primus/mfdata/"
    target_dir = "/global/data/scr/chh327/primus/data/target/"
    
    mfdata_fname = "mfdata_all_supergrid01_sdss_lit.fits.gz"
    target_noedgecut_fname = "target_cylr2h50_thresh75_nbin5_all_sdss_lit.fits"

    mfdata = mrdfits(mfdata_dir+mfdata_fname,1)
    target_noedgecut = mrdfits(target_dir+target_noedgecut_fname,1)
   
;    plot, mfdata.z, mfdata.mass, psym=3 
;    oplot, target_noedgecut.redshift, target_noedgecut.mass, psym=1
;    print, min(target_noedgecut.mass), max(target_noedgecut.mass)
    
    mfdata_smf = im_mf_vmax(mfdata.mass,mfdata.weight/mfdata.vmax_evol,masslimit=mfdata.masslimit)
    target_smf = im_mf_vmax(target_noedgecut.mass,target_noedgecut.weight/target_noedgecut.vmax,masslimit=target_noedgecut.masslimit)
    target_incorrect_smf = im_mf_vmax(target_noedgecut.mass,target_noedgecut.weight/target_noedgecut.vmax_evol,masslimit=target_noedgecut.masslimit)
    target_noedgecut_smf = im_mf_vmax(target_noedgecut.mass,target_noedgecut.weight/target_noedgecut.vmaxavail,masslimit=target_noedgecut.masslimit)
    mwrfits,mfdata_smf,"/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/mfdata_smf.fits",/create
    mwrfits,target_smf,"/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/target_smf.fits",/create
    mwrfits,target_incorrect_smf,"/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/target_incorrect_smf.fits",/create
    mwrfits,target_noedgecut_smf,"/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/target_noedgecut_smf.fits",/create
end
