import shutil
from ipy             import *
from FPI_util_funcs  import *
from OpPt_util_funcs import *

unit        = 'FM203'
head0       = 'H0'
head1       = 'H1'
start_orb   = 590
end_orb     = 1000
obs, bay    = serial_num_to_sc(unit)

Thrsh_V_str = '%s_%s_ThrshV'       % (obs, bay)
Counts0_str = '%s_%s_counts_Head0' % (obs, bay)
Counts1_str = '%s_%s_counts_Head1' % (obs, bay)

base_dir = 'x:/data/ftp/mms2/fpi/cal/OpPtCal/'
my_list1 = find_oppt_CDFs(base_dir,start_orb,end_orb,unit)

Energy_list = ['Energy'+str(i) for i in range(8)]
for m in my_list1:
    skippy  = pycdf.CDF(m)
    mode    = skippy.attrs[22][0]
    ThrshV  = np.asarray(skippy[Thrsh_V_str])
    Counts0 = np.asarray(skippy[Counts0_str])
    Counts1 = np.asarray(skippy[Counts1_str])
    skippy.close()
    SweepV  = ThrshV[0:16]
    C0      = create_OpPt_Dict(Counts0)
    C1      = create_OpPt_Dict(Counts1)
    orbit   = get_orb_number(m)
    for e in Energy_list:
        title = 'Pixel Sweep Plot for %s%s on orbit %s for energy %s at %s mode' % (unit, head0, orbit, e, mode)
        fname = 'c:/users/cschiff/python/FitAgnostic/%s%s_Sweep_%s_%s_%s.png' % (unit, head0, orbit, e, mode)
        stat_plot_energy(SweepV,C0,e,title,fname)
        title = 'Pixel Sweep Plot for %s%s on orbit %s for energy %s at %s mode' % (unit, head1, orbit, e, mode)
        fname = 'c:/users/cschiff/python/FitAgnostic/%s%s_Sweep_%s_%s_%s.png' % (unit, head1, orbit, e, mode)
        stat_plot_energy(SweepV,C1,e,title,fname)