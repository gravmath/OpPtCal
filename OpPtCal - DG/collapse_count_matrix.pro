FUNCTION collapse_count_matrix,specie,curser,curswp

if specie eq 'DES' then begin
  counts = total(total(total(curser[curswp].counts, 3), 3), 3)
endif
  
  
if specie eq 'DIS' then begin
  counts_sum_sweep = total(total(curser[curswp].counts, 3), 3)
  counts1 = counts_sum_sweep[*,*,0] + counts_sum_sweep[*,*,1]
  counts2 = counts_sum_sweep[*,*,1] + counts_sum_sweep[*,*,7]
  counts3 = counts_sum_sweep[*,*,1] + counts_sum_sweep[*,*,6]
  counts4 = counts_sum_sweep[*,*,2] + counts_sum_sweep[*,*,6]
  counts5 = counts_sum_sweep[*,*,2] + counts_sum_sweep[*,*,5]
  counts6 = counts_sum_sweep[*,*,3] + counts_sum_sweep[*,*,4]
  counts_array = [total(counts1),total(counts2),total(counts3),total(counts4),total(counts5),total(counts6)]
  result = where(counts_array eq max(counts_array))
  
    if result eq 0 then counts = counts1
    if result eq 1 then counts = counts2
    if result eq 2 then counts = counts3
    if result eq 3 then counts = counts4
    if result eq 4 then counts = counts5
    if result eq 5 then counts = counts6
  
;  if result eq 0 then print, 'Max Counts in Energies: 500-3000 eV' & counts = counts1
;  if result eq 1 then print, 'Max Counts in Energies: 1500-3000 eV' & counts = counts2
;  if result eq 2 then print, 'Max Counts in Energies: 3000-5000 eV' & counts = counts3
;  if result eq 3 then print, 'Max Counts in Energies: 5000-7000 eV' & counts = counts4
;  if result eq 4 then print, 'Max Counts in Energies: 7000-10000 eV' & counts = counts5
;  if result eq 5 then print, 'Max Counts in Energies: 15000-20000 eV' & counts = counts6
endif

return, counts
END