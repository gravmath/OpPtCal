; auxiliary functions for thrswp_fit.pro

FUNCTION phd_cdf, tr,Gamma,Q
  ; INPUT
  ;   tr - preamplifier threshold
  ;   Gamma - MCP pulse height distribution (PHD) shape
  ;   Q - MCP PHD scale 
  ; OUTPUT
  ;   cumulative distribution function of MCP PHD at threshold tr
  
  IF gamma lt 0 then gamma = 0.
  f = 1d0 - IGAMMA(Gamma,tr/Q)
  RETURN, f
END

FUNCTION getClippedCounts_AllPX, tr,lim0CRs,Gammas,Qs,AAXT,px
  ; INPUT
  ;   tr - preamplifier threshold
  ;   lim0CRs[px0...px15] - limit t->0 countrate for each pixel (convolved MCP efficiency w/ sampled PSD)
  ;   Gammas[px0...px15] - Gamma parameter for gamma distribution for each pixel
  ;   Qs[px0...px15] - Q parameter for gamma distribution for each pixel
  ;   AAXT[2] -Anode-Anode X-talk 1 pixel and 2 pixels away
  ;   px - pixel to accumulate counts  
  ; OUTPUT
  ;   Measured counts in a pixel Px at threshold tr
 
  ; Start with no anode-anode crosstalk losses - increase loss for adjacent pixels
  AAXT0 = 1d0

  ; Contribution of PX-2 to PX
  IF px GT 1d0 THEN BEGIN
    Niim2 = lim0CRs[px-2]*phd_cdf(tr/AAXT[1],Gammas[px],Qs[px])
    AAXT0 = AAXT0 - AAXT[1]
  ENDIF ELSE BEGIN 
    Niim2 = 0d0
  ENDELSE
  
  ; Contribution of PX-1 to PX
  IF px GT 0d0 THEN BEGIN 
    Niim1 = lim0CRs[px-1]*phd_cdf(tr/AAXT[0],Gammas[px],Qs[px])
    AAXT0 = AAXT0 - AAXT[0]
  ENDIF ELSE BEGIN 
    Niim1 = 0d0
  ENDELSE 
  
  ; Contribution of PX+1 to PX
  IF px LT 15d0 THEN BEGIN 
    Niip1 = lim0CRs[px+1]*phd_cdf(tr/AAXT[0],Gammas[px],Qs[px])
    AAXT0 = AAXT0 - AAXT[0]
  ENDIF ELSE BEGIN 
    Niip1 = 0d0
  ENDELSE 
  
  ; Contribution of PX+2 to PX
  IF px LT 14d0 THEN BEGIN 
    Niip2 = lim0CRs[px+2]*phd_cdf(tr/AAXT[1],Gammas[px],Qs[px]) 
    AAXT0 = AAXT0 - AAXT[1]
  ENDIF ELSE BEGIN 
    Niip2 = 0d0
  ENDELSE

  ; Contribution of PX to PX
  Nii = lim0CRs[px]*phd_cdf(tr/aaxt0,Gammas[px],Qs[px])
  
  ; Return log10 of total counts in PX
  RETURN, ALOG10(Niim2+Niim1+Nii+Niip1+Niip2)
END



FUNCTION getClippedCounts_SingPX, tr,lim0CR,Gamma,Q,AAXT,scale_factor
  ; INPUT
  ;   tr - preamplifier threshold
  ;   lim0CR - limit t->0 countrate for  pixel (convolved MCP efficiency w/ sampled PSD)
  ;   Gamma - Gamma parameter for gamma distribution for pixel
  ;   Q - Q parameter for gamma distribution for pixel
  ;   AAXT[3] -Anode-Anode X-talk 1 pixel and 2 pixels away and PIXEL to use
  ; OUTPUT
  ;   Measured counts in a pixel Px at threshold tr
 
  ; Start with no anode-anode crosstalk losses - increase loss for adjacent pixels
  AAXT0 = 1d0
  px = AAXT[2]

  ; Contribution of PX-2 to PX
  IF px GT 1d0 THEN BEGIN
    Niim2 = lim0CR*phd_cdf(tr/AAXT[1],Gamma,Q)
    AAXT0 = AAXT0 - AAXT[1]
  ENDIF ELSE BEGIN 
    Niim2 = 0d0
  ENDELSE
  
  ; Contribution of PX-1 to PX
  IF px GT 0d0 THEN BEGIN 
    Niim1 = lim0CR*phd_cdf(tr/AAXT[0],Gamma,Q)
    AAXT0 = AAXT0 - AAXT[0]
  ENDIF ELSE BEGIN 
    Niim1 = 0d0
  ENDELSE 
  
  ; Contribution of PX+1 to PX
  IF px LT 15d0 THEN BEGIN 
    Niip1 = lim0CR*phd_cdf(tr/AAXT[0],Gamma,Q)
    AAXT0 = AAXT0 - AAXT[0]
  ENDIF ELSE BEGIN 
    Niip1 = 0d0
  ENDELSE 
  
  ; Contribution of PX+2 to PX
  IF px LT 14d0 THEN BEGIN 
    Niip2 = lim0CR*phd_cdf(tr/AAXT[1],Gamma,Q) 
    AAXT0 = AAXT0 - AAXT[1]
  ENDIF ELSE BEGIN 
    Niip2 = 0d0
  ENDELSE

  ; Contribution of PX to PX
  Nii = lim0CR*phd_cdf(tr/aaxt0,Gamma,Q)

  ; Return log10 of total counts in PX
  RETURN, ALOG10(Niim2*scale_factor+Niim1*scale_factor+Nii+Niip1*scale_factor+Niip2*scale_factor)
END



FUNCTION tswp_AllPX,X,P,AAXT=AAXT
; INPUT: 
;   X[16xT] - for each of 16 pixels have T thresholds [tr0px0,...,trTpx0,tr0px1,...,trTpx15]
;   P[48]
;     0-15:     lim0CR[px0...px15] - limit t->0 countrate for each pixel (convolved MCP efficiency w/ sampled PSD)
;     16-31:    Gamma[px0...px15] - Gamma parameter for gamma distribution for each pixel
;     32-47:    Q[px0...px15] - Q parameter for gamma distribution for each pixel
;   AAXT[2]
;     Anode-Anode X-talk 1 pixel and 2 pixels away  
; OUTPUT: 
;     Y[16xT] - log10 of counts for each pixels at each threshold
 
  lim0CRs = P[0:15]
  Gammas = P[16:31]
  Qs =  P[32:47]
  
  T = size(X,/N_ELEMENTS)/16 
  Y = DBLARR(size(X,/N_ELEMENTS))
  
  IF (n_elements(WHERE(P LE 0,/null)) EQ 0) THEN BEGIN
    FOR px = 0,15 DO BEGIN   
      trs = X[(px*T):(px*T+T-1)]
      Y[(px*T):(px*T+T-1)] = getClippedCounts_AllPX(trs,lim0CRs,Gammas,Qs,AAXT,px)
    ENDFOR
  ENDIF
  
  Y[WHERE(X EQ 0)] = 0d0
  RETURN, Y
END

 
FUNCTION tswp_SingPX,X,P,AAXT=AAXT
; INPUT: 
;   X[T] - for each pixel have T thresholds [tr0,...,trT]
;   P[4] - lim0CR, Gamma, Q, scale_factor
;   AAXT[3]- Anode-Anode X-talk 1 pixel [0] and 2 pixels away [1] and [2] current pixel to use  
; OUTPUT: 
;     Y[T] - log10 of counts at each threshold
 
  lim0CR = P[0]
  Gamma = P[1]
  Q =  P[2]
  scale_factor = P[3]
  
  ;T = size(X,/N_ELEMENTS)
  Y = DBLARR(size(X,/N_ELEMENTS))
  
  IF (n_elements(WHERE(P LE 0,/null)) EQ 0 and scale_factor lt 10) THEN BEGIN
      Y = getClippedCounts_SingPX(X,lim0CR,Gamma,Q,AAXT,scale_factor)
  ENDIF
  
  Y[WHERE(X EQ 0)] = 0d0
  Y[where(FINITE(y) eq 0)] = 0.
  RETURN, Y
END


