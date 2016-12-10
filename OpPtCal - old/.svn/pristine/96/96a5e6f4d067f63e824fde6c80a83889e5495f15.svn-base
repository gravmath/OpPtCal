; Run OpPntCal proc.
;
; Last modified: MPH, 13 Jan, 2015.
;
PRO runOpPntCal_on_folder ; No params.

;FPIprocToolsDir = 'C:\fpi\FPIdataProcTools\'
;defSysV, '!FPIthresholdConvert_PATH', $ ; Absolute path to:
;           FPIprocToolsDir + $          ;  FPI threshold-
;          'FPIthresholdConvert\', 1     ;  conversion data (Read-Only).
;   print, !FPIthresholdConvert_PATH
;
;thrSwp_func ; Force compilation.
;defSysV, '!FPIGroundCal_PATH', $        ; Absolute path to:
;           FPIprocToolsDir + $          ;  Cal-analysis tools (R/O).
;          'CalAnalysis\GroundCal\', 1
;   print, !FPIGroundCal_PATH
;
;device, Decomposed=0 ; Set to plot w/ 8-bit colors.
;loadCT, 39, /Silent ; Load Rainbow+White ColorTable.
;
;;root = 'C:\Users\fpi_egse\Desktop\FPI_MRT17_review\MRT17_HV_TurnOn\' + $
;;       'RevE0000_00\FPI005_Q3_Ramp_and_Checkout_00\raw\'

argv = COMMAND_LINE_ARGS(COUNT=argc)

root = argv[0]

processOpPtCal, root ; Run proc.

end
