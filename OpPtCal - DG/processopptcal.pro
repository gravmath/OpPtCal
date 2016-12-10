pro processOpPtCal, session, stub=stub, csv=csv

if n_elements(stub) eq 0 then stub=0
if n_elements(csv) eq 0 then csv=0

fit_allPX = 'no'

print, 'Checking for OpPtCal in: ' + session
outdir = strmid(session, 0, strlen(session) - 4) + 'analyzed' + path_sep()

data = FPIcalDataPkg_OperPnt_fetch(session, stub=stub)
;restore, '~/opptcal.sav'

colors=['blue', 'red', 'gold', 'pur5', 'orange']

if data eq !NULL then begin
  cgDisplay, 1500, 900, /pixmap
  cgtext, .1, .4, 'DATA FILE COULD NOT BE READ!!', charsize=4
  specie='XXX'
  ser=000
  ttl='OpPtCal ' + specie + string(ser, format='(I03)')
  write_png, strrep(outdir + ttl + stripslashes(session) + '.png', ' ', '_'), tvrd(/true)
  
endif

save, data, filename=outdir+'IDLData.sav'

;pos = cgLayout([4, 4], OXMargin=[2, 2], OYMargin=[50, 3], YGap=0, xgap=0)
;pos2= cglayout([2,2],  OXMargin=[6, 6], OYMargin=[5, 4], YGap=3, xgap=12)

pos = alayout([4,4], ox=[0.02, 0.02], oy=[0.50, 0.03], xg=0, yg=0)
pos2= alayout([2,1], ox=[0, 0], oy=[.08, .53],  xg=.1, yg=0)

x=data.serial

w=0


;foreach ser, x[uniq(x[sort(x)])] do begin ;i dont know why IDL requires the data to be sorted in order for uniq to work.....
xsort = x[sort(x)]  ; changed by DJG 04/13/2015 - not sure why this was required but above line was only giving us DIS data...
foreach ser,xsort[uniq(xsort)] do begin  
  cgdisplay, 1500, 900, wid=w, xpos=w*1200, ypos=0, /pixmap
  w=w+1

  if ser lt 100 then specie='DIS' else specie='DES'
  erase, 'ffffff'x
  print, 'Fitting data for ' + string(ser, format='(I03)')
  
  gndh0 = read_specfit(ser, 0)
  gndh1 = read_specfit(ser, 1)
  
  gndgain=[[gndh0.gains], [gndh1.gains]]
  gndetag=[[IGAMMA(gndh0.gammas,4.0e5/gndh0.Qs)], [IGAMMA(gndh1.gammas,4.0e5/gndh1.Qs)]]
  
  ;; sum over energy/iteration
  curser = data[where(data.serial eq ser)]
  
  nels = n_elements(curser)
  cgLoadCT, 34, ncolors=nels
  cgtext, .6 + 0.040 * nels, .01, 'CurrentCal', color='grn5', /normal
  cgtext, .1 + 0.040 * nels, .01, 'CurrentCal', color='grn5', /normal
    
  for curswp = 0, n_elements(curser) - 1 do begin

    cgtext, .1 + 0.040 * curswp, .01, string(curser[curswp].MCPvlt[0], format='(I05)') + 'V', color=colors[curswp], /normal
    cgtext, .6 + 0.040 * curswp, .01, string(curser[curswp].MCPvlt[1], format='(I05)') + 'V', color=colors[curswp], /normal
    ;print, string(curser[curswp].MCPvlt[0]) + 'V'
    for head = 0,1 do begin
      print, 'Processing ' + string(ser, format='(I03)') + ', H' + string(head, format='(I1)') + ', MCPB=' + string(curser[curswp].mcpvlt)
      ;counts = total(total(total(curser[curswp].counts, 3), 3), 3)
      counts=reform((total(total(curser[curswp].counts, 4), 4))[*,*,head])
      
      aaxt   = curser[curswp].anodXT
      thresh = curser[curswp].thresh[*,head]
      
        ;;;;; just for testing;;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;
        
          if stub then thresh = thresh * ((curser[curswp].mcpvlt / max(curser.mcpvlt))[head])^5 * 1.6
        
        ;;;;;remove later;;;;;;;;; 
        ;;;;;;;;;;;;;;;;;;;;;;;;;;
      tfit = thrswp_fit(TRANSPOSE(counts), thresh, aaxt,specie,fit_allPX)
      etaGs = IGAMMA(tfit.gammas,4.0e5/tfit.Qs)
      
      trs_model  = DINDGEN(100)/99.0*5 + 3 ; generate thresholds between 1e3 and 1e8 for modeled curve
      trs_model = 10^trs_model 
  
      ;;;plot th sweeps
      
      for pix = 0, 15 do begin
        
        ;account for bad fit
        if tfit.lim0crs[pix] eq 0 then tfit.lim0crs[pix] = max(counts[*,pix])
        
        if fit_allPX eq 'no' then begin
          N_model_px = getClippedCounts_SingPX(trs_model, tfit.lim0crs[pix], tfit.gammas[pix], tfit.qs[pix], [aaxt[0],aaxt[1],pix],tfit.scale_factor[pix])
          ;print,'Scale Factor: '+string(tfit.scale_factor[pix])
        endif else begin
           N_model_px = getClippedCounts_allPX(trs_model, tfit.lim0crs, tfit.gammas, tfit.qs, aaxt,pix)  ; obtain modeled curve from fit parameters
        endelse
        
        N_model_px = 10^N_model_px
         
        cgplot, thresh, counts[*,pix]/tfit.lim0crs[pix], color=colors[curswp], psym=2, /ylog, /xlog, $
                 position = pos[*,pix]/[2.0, 1, 2.0, 1] + [0.5, 0, 0.5, 0] * head, /noerase, $
                 xrange=[5e4, 5e7], yrange=[.1, 4], XTICKFORMAT="(A1)", YTICKFORMAT="(A1)"
                 
        cgplot, trs_model, N_model_px/tfit.lim0crs[pix], color=colors[curswp], /overplot, thick=2
        cgplot, [tfit.gains[pix], tfit.gains[pix]], [.1, 10], color=colors[curswp], /overplot, thick=2
      
        cgplot, [2e6, 2e6], [.1, 10], color='black', thick=2, /overplot, linestyle=2
        cgtext, pos[0,pix]/2.0+.01+0.5*head, pos[1, pix]+.01, string(pix, format='(I02)'), /normal  
      endfor
      
      ;; plot gains
      
      err = [tfit.gains + tfit.gains_err, reverse(tfit.gains - tfit.gains_err), tfit.gains[0] + tfit.gains_err[0]]
      ind = [indgen(16),  reverse(indgen(16)), 0]
      
      cgplot, [0,15], [2e6, 2e6], color='black', /nodata, yrange=[5e5, 5e7], /noerase, position=pos2[*,head], /normal, /ylog, ystyle=8, ytitle='Gain', xtitle='Pixel'
      
      ;cgcolorfill, ind, err, color=colors[curswp]
      
      cgplot, indgen(16), tfit.gains, color=colors[curswp], thick=3, yrange=[5e5, 5e7], /overplot, position=pos2[*,head], /normal, /ylog
         
      cgplot, [0,15], [2e6, 2e6], color='black', thick=2, /overplot, linestyle=2
      
      ;put in ground cal
      cgplot, indgen(16), gndgain[*,head], color='grn5', thick=4, /overplot

      ;signal loss
      cgaxis, yaxis=1, yrange=[0, 25],  /save, ytitle='(---) Signal Loss (%)', ylog=0, xtitle='Pixel'
      cgplot, indgen(16), etaGs*100, color=colors[curswp], thick=2, linestyle=2, /overplot
      cgplot, indgen(16), gndetag[*,head], color='grn5', thick=2, linestyle=2, /overplot
      
      ;;;;;output csvs;;;;;
      ;
      
      if csv then begin
      
        ;gain
      
        csvstr = strjoin([string(curser[curswp].strtjd, format = '(F9.1)'), string(curser[curswp].MCPVLT[head], format='(I)'), string(tfit.gains), string(tfit.gains_err)], ',')
        csvfil = specie + string(ser, format='(I03)') + 'H' + string(head, format='(I1)') + 'MCP' + string(curser[curswp].MCPVLT[head], format='(I04)') + '_Gain.csv
        dirname=!CalCSVDir
        openw, lun, dirname + path_sep() + csvfil, /get_lun, /append
        printf, lun, csvstr
        close, lun
        free_lun, lun
        
        ;sig loss
        
        csvstr = strjoin([string(curser[curswp].strtjd, format = '(F9.1)'), string(curser[curswp].MCPVLT[head], format='(I)'), string(etaGs)], ',')
        csvfil = specie + string(ser, format='(I03)') + 'H' + string(head, format='(I1)') + 'MCP' + string(curser[curswp].MCPVLT[head], format='(I04)') + '_SignalLoss.csv
        dirname=!CalCSVDir
        openw, lun, dirname + path_sep() + csvfil, /get_lun, /append
        printf, lun, csvstr
        close, lun
        free_lun, lun
        
      endif
    endfor
    
    ;stop
  endfor
  
          
  cgtext, .05, .01, 'Head 0', /normal
  cgtext, .55, .01, 'Head 1', /normal
  
  cgtext, 0.5025, pos[1,12], 'Normalized Count Rate', orientation=90, /normal, charsize=1
  cgtext, .51, .488, 'Threshold', /normal, charsize=1
  ;cgcolorfill, 

  ttl='OpPtCal ' + specie + string(ser, format='(I03)')
  cgtext, 0.05, 0.975, ttl, /normal, charsize=2
  cgtext, 0.3, 0.975, session, charsize=1, /normal
      
  print, 'writing - ' +  strrep(outdir + ttl + stripslashes(strrep(session, ':', '_')) + '.png', ' ', '_')
  write_png, strrep(outdir + ttl + stripslashes(strrep(strmid(session, strpos(session, 'FPI')), ':', '_')) + '.png', ' ', '_'), tvrd(/true)
  
endforeach

end

