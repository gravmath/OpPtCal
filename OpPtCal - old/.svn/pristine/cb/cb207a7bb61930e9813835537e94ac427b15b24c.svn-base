;+
;  :Description:
;     This function performs a non-linear least squares fit on a threshold sweep collected for a set of pixels
;     
;  :Categories:
;     FPI Calibration Data Proc
;
;  :Authors:
;     Daniel J. Gershman
;
;  :History:
;     08/14/2014 -- Daniel J. Gershman, NASA-GSFC:  initial version
;
;  :Params:
;
;     N : in, type=dblarr(16, T), counts from threshold sweep for P pixels at T thresholds
;
;     thresholds: in, type=dblarr(T),  array of T thresholds for sweep
;     
;     aaxt:   in, type=dblarr(2),      anode-anode crosstalk for adjacent pixels - index 0 and 1 correspond to pixels distance +-1 and +-2, respectively
;
;     tfit : out, 
;        struct = {lim0crs, dblarr(16),          $     ;   number of counts for tr->0 for each pixel 0...15
;                  lim0crs_err, dblarr(16),      $     ;   1-sigma uncertainty in lim0cr for each pixel 
;                  gammas, dblarr(16),           $     ;   gamma parameter for MCP pulse height distribution (PHD) for each pixel
;                  gammas_err, dblarr(16),       $     ;   1-sigma uncertainty in gamma parameter for each pixel
;                  Qs, dblarr(16),               $     ;   Q parameter for MCP PHD for each pixel
;                  Qs_err, dblarr(16),           $     ;   1-sigma uncertainty in Q parameter for each pixel
;                  Gains, dblarr(16),            $     ;   median gain for MCP PHD for each pixel
;                  Gains_err, dblarr(16),        $     ;   1-sigma uncertainty in median gain for each pixel
;                  V, dblarr(48,48)              $     ;   covariance matrix (symmetric) of least-squares fit  
;                     index 0-15:     lim0crs[px0...px15] 
;                     index 16-31:    gammas[px0...px15] 
;                     index 32-47:    Qs[px0...px15]   
;                 R2, dblarr(1),                       ; value of the summed squared residuals
;                 }
;                 
;-


FUNCTION thrswp_fit, N, thresholds, aaxt,specie,fit_allPX

      ; Initialize fit structure
    tfit = {  lim0crs: dblarr(16),          $     ;   number of counts for tr->0 for each pixel 0...15
              lim0crs_err: dblarr(16),      $     ;   1-sigma uncertainty in lim0cr for each pixel 
              gammas: dblarr(16),           $     ;   gamma parameter for MCP pulse height distribution (PHD) for each pixel
              gammas_err: dblarr(16),       $     ;   1-sigma uncertainty in gamma parameter for each pixel
              Qs: dblarr(16),               $     ;   Q parameter for MCP PHD for each pixel
              Qs_err: dblarr(16),           $     ;   1-sigma uncertainty in Q parameter for each pixel
              Gains: dblarr(16),            $     ;   median gain for MCP PHD for each pixel
              Gains_err: dblarr(16),        $     ;   1-sigma uncertainty in median gain for each pixel
              V: dblarr(48,3),             $     ;   covariance matrix (symmetric) of least-squares fit  
              R2: dblarr(1),                $     ;   value of the summed squared fit residual
             fit_flag: dblarr(16),          $     ;   flag to report wheter high, mid, and low gain
             scale_factor: dblarr(16)}            ;   scale factor for the single pixel fits
    
    CATCH, Error_status
    
    if specie eq 'DES' then aaxt = [0.022,7e-3,0] 
    if specie eq 'DIS' then aaxt = [2.5e-3,1e-6,0] 
    ;;print, 'AAXT = ' +string(aaxt)
    
    IF Error_status NE 0 THEN BEGIN
      PRINT, 'Error index: ', Error_status
      PRINT, 'Error message: ', !ERROR_STATE.MSG
      ; Handle the error by extending A:
      return, tfit
      CATCH, /CANCEL
    ENDIF
   
  ; Format data for fit
    
    
    T = N_ELEMENTS(thresholds)
          
    X = DBLARR(16*T)
    Y = DBLARR(16*T)

    FOR px = 0,15 DO BEGIN   
      X[(px*T):(px*T+T-1)] = thresholds
      Y[(px*T):(px*T+T-1)] = (ALOG10(N[px,*]))
    ENDFOR
 
  ; set X and Y = 0 for invalid data
    X[WHERE(FINITE(Y) EQ 0)] = 0d0
    Y[WHERE(FINITE(Y) EQ 0)] = 0d0
   
  ; Initial guesses for parameters
    P_allPX = DBLARR(48)
    P_singPX = DBLARR(3)
        
    iiT = WHERE(thresholds GE 4.0e5)
    
    

FOR px = 0,15 DO BEGIN
      pxx = floor((px)/4.)
      pxy = string(px)
      Nd00 = N[px,*]
      
      ;;; initial guesses
      gamma = 1.5
      g_sig = 5e6
      g_xt = aaxt[0]*5e6
      lim0cr = 0
       
      
  if total(Nd00) gt 0 then begin
         
    p1 = Nd00[0:4]
    p2 = Nd00[5:9]
    p3 = Nd00[10:-1]
    mms = [mean(p1),mean(p2),mean(p3)]
    mss = [stddev(p1)/mean(p1),stddev(p2)/mean(p2),stddev(p3)/mean(p3)]
    ii = (where(mss eq min(mss),/null))

    if ii eq !NULL then ii = 0 else ii = ii[0]
    
    tfit.fit_flag[px]=3

    ;;print, mms
    ;;print, mss
    ; compare variations in first, middle, and last third of data
    
    if ii eq 0 and mss[0]/(mss[2]>1e-9) gt 0.5 then begin
        ; if both the first and last third are each not varying strongly, compare their ratio to distinguish between low and high gain states
          lim0cr = mms[2]
          tfit.fit_flag[px] = 2 
    endif else begin
      
      if ii eq 2 and mss[2] eq 0 then begin
        lim0cr = mms[0] 
        tfit.fit_flag[px]=1 
        endif else begin 
          lim0cr = mms[ii]
          if ii eq 0 then tfit.fit_flag[px]=2
          if ii eq 1 then tfit.fit_flag[px]=1
          if ii eq 2 then tfit.fit_flag[px]=0
        endelse  
    endelse

    if lim0cr eq !NULL then lim0cr = 0
    
    ;help, lim0cr
    
        g_xt = where(Nd00 gt lim0cr[0]*2.,/null)
        g_sig = where(Nd00 lt lim0cr[0]*0.5,/null)
         
        if_dummy = 0
        if n_elements(g_sig) gt 0 and n_elements(g_xt) gt 0 and if_dummy eq 0 then begin
          gamma = 1.5
          g_sig = thresholds[g_sig[0]]
          g_xt = g_sig*aaxt[0]
          if_dummy = 1
        endif
        
        if n_elements(g_sig) gt 0 and n_elements(g_xt) eq 0 and if_dummy eq 0 then begin
          gamma = 1.
          g_sig = thresholds[g_sig[0]]
          g_xt = g_sig*aaxt[0]
          if_dummy = 1
        endif
        
        if n_elements(g_sig) eq 0 and n_elements(g_xt) gt 0 and if_dummy eq 0 then begin
          gamma = 2.
          g_xt = thresholds[g_xt[-1]]
          g_sig = g_xt/aaxt[0]
          if_dummy = 1
        endif

        if if_dummy eq 0 then begin
          gamma = 1.5
          g_sig = 5e6
          g_xt = aaxt[0]*5e6
        endif
        
         

   endif
      
      Q = g_sig/(gamma*(3.*gamma-0.8)/(3.*gamma+0.2)) 
    
      
      ;iiTN = WHERE(TS_DIFF(TRANSPOSE(N[px,iiT]),1) EQ MAX(TS_DIFF(TRANSPOSE(N[px,iiT]),1)))
      ;Qguess = thresholds[iiT[iiTN[0]]]
      P_allPX[px] = lim0cr[0]
      P_allPX[16+px] = gamma[0]
      P_allPX[32+px] = Q[0]      ; guess for Q
      
      P_singPX = [lim0cr[0],gamma[0],Q[0],1.]
      
      ;help, P_singPX
      
      if fit_allPX eq 'no' then begin

      PARAM_INFO = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
      res = 0d0
      V = MAKE_ARRAY(4,4,/DOUBLE,VALUE=0)      
      AAXT[2] = px
      
      cX = thresholds
      cY = (ALOG10(N[px,*]))
      cX[WHERE(FINITE(cY) EQ 0)] = 0d0
      cY[WHERE(FINITE(cY) EQ 0)] = 0d0
      
      
      fit_params = MPFITFUN('tswp_SingPX',cX,cY,DBLARR(N_ELEMENTS(cY)),WEIGHTS=cY,P_singPX,FUNCTARGS={AAXT:AAXT},COV=V,/QUIET,/DOUBLE,PARINFO=PARAM_INFO,BESTNORM=res)
      DOF = N_ELEMENTS(cY) - N_ELEMENTS(fit_params) ; deg of freedom
      V = V/DOF

      tfit.V[px*3+0,0] = V[0,0]
      tfit.V[px*3+0,1] = V[0,1]
      tfit.V[px*3+0,2] =  V[0,2]
      tfit.V[px*3+1,0] =  V[1,0]
      tfit.V[px*3+1,1] =  V[1,1]
      tfit.V[px*3+1,2] =  V[1,2]
      tfit.V[px*3+2,0] =  V[2,0]
      tfit.V[px*3+2,1] =  V[2,1]
      tfit.V[px*3+2,2] =  V[2,2]

      
      ;tfit.V[px,px] = V[0]
      ;tfit.V[16+px,16+px] = V[1]
      ;tfit.V[32+px,32+px] = V[2]

      tfit.R2 = tfit.R2 + res
      
      tfit.lim0crs[px] = fit_params[0]
      tfit.gammas[px] = fit_params[1]
      tfit.Qs[px] = fit_params[2]
      tfit.scale_factor[px] = fit_params[3]

      endif

ENDFOR
    
    

       
  ; Fit data - weight by the log10 counts (i.e., Y)

;   if fit_allPX eq 'yes' then begin
;    ;print,'Fitting all pixels!'
;     PARAM_INFO = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},48)
 ;    res = 0d0
;     V = MAKE_ARRAY(48,48,/DOUBLE,VALUE=0)


;    fit_params = MPFITFUN('tswp_AllPX',X,Y,DBLARR(N_ELEMENTS(Y)),WEIGHTS=Y,P_allPX,FUNCTARGS={AAXT:AAXT},COV=V,/QUIET,/DOUBLE,PARINFO=PARAM_INFO,BESTNORM=res)
;    DOF = N_ELEMENTS(Y) - N_ELEMENTS(fit_params) ; deg of freedom
;    V = V/DOF
    ;tfit.V = V
;    tfit.R2 = res

 ;   tfit.lim0crs = fit_params[0:15]
 ;   tfit.gammas = fit_params[16:31]
 ;   tfit.Qs = fit_params[32:47]


 ;  endif


  ; Construct covariance matrix (normalized by DOF)




  ; extract out fit parameters and MCP gain
    
    for pixel=0,15 do begin
      if tfit.gammas[pixel] gt 1 then begin
        ; for gamma > 1, can use approximate formula for median
        tfit.Gains[pixel]=tfit.Qs[pixel]*tfit.gammas[pixel]*(3.0*tfit.gammas[pixel]-0.8)/(3.0*tfit.gammas[pixel]+0.2)
      endif else begin
        ; for gamma <= 1, use exponential approximation for median
        tfit.Gains[pixel]=tfit.Qs[pixel]*alog(2.)
      endelse 
    endfor
    
    
    ;tfit.Gains = tfit.Qs*tfit.gammas*(3.0*tfit.gammas-0.8)/(3.0*tfit.gammas+0.2)
    
  ; estimate 1sig fit errors from diagonal of normalized covariance matrix

    ;errs = SQRT(DIAG_MATRIX(tfit.V))
    for pixel=0,15 do begin
      tfit.lim0crs_err[pixel] = SQRT(tfit.V[3*pixel+0,0])
      tfit.gammas_err[pixel] = SQRT(tfit.V[3*pixel+1,1])
      tfit.Qs_err[pixel] = SQRT(tfit.V[3*pixel+2,2])
    endfor

    ;tfit.lim0crs_err = errs[0:15]
    ;tfit.gammas_err = errs[16:31]
    ;tfit.Qs_err = errs[32:47]
    
;;print, tfit.lim0crs_err/tfit.lim0crs
;;print, tfit.gammas_err/tfit.gammas
;;print, tfit.Qs_err/tfit.Qs

  ; estimate median gain 1sig errors - use generalized uncertainty formula to estimate (linear) errors
    Vij = DBLARR(16)
    FOR ii=0,15 DO BEGIN
      ;Vij[ii] = tfit.V[16+ii,32+ii]
      Vij[ii] = tfit.V[3*ii+1,2]  ; covariance between gamma and Q
    ENDFOR

;;print, Vij

    dGdGamma = tfit.Qs*(((3.0*tfit.gammas-0.8)*(3.0*tfit.gammas+0.2)+3)/(3.0*tfit.gammas+0.2)^2)
    dGdQ = tfit.gammas*(3.0*tfit.gammas-0.8)/(3.0*tfit.gammas+0.2)
    tfit.Gains_err = SQRT( (dGdQ*tfit.Qs_err)^2 + (dGdGamma*tfit.gammas_err)^2 + 2*dGdQ*dGdGamma*Vij)
;;print, tfit.Gains_err/tfit.Gains
    return, tfit
END