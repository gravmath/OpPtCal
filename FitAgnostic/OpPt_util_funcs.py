from   ipy import *
from   FPI_util_funcs import *
import matplotlib.pyplot as plt

def period_fold(array,stride):
    num_rows   = len(array)
    num_slices = num_rows/stride 
    pf_array   = np.zeros(stride)
    for i in range(stride):
        pf_array[i] = np.sum(array[i::stride])
    
    return pf_array
    
def create_OpPt_Dict(Counts_array):
    head_dict  = {}
    for pixel in range(0,16):
        if pixel < 10:
            pixel_string = '0'+str(pixel)
        else:
            pixel_string = str(pixel)
        Pixel_label = 'Pixel'+pixel_string 
        head_dict[Pixel_label] = {'Energy0':[],'Energy1':[],'Energy2':[],'Energy3':[],\
                                  'Energy4':[],'Energy5':[],'Energy6':[],'Energy7':[]}

        for iteration in range(0,8):
            for sweep in range(0,8):
                energy_index  = 'Energy'+str(sweep)
                start_slice   = 128*iteration + 16*sweep
                stop_slice    = 128*iteration + 16*(sweep+1)
                head_dict[Pixel_label][energy_index].append(Counts_array[start_slice:stop_slice,pixel])
    
    return head_dict
    
def stat_plot_energy(ThreshV,head_dict,energy,title_string,plot_name):
    fig, axes  = plt.subplots(4,4,figsize=(16,16))
    colors = ['bisque','chartreuse','aliceblue','darkmagenta','aquamarine','darkturquoise','crimson','fuchsia']
    colors = ['blue','red','green','cyan','magenta','darkorchid','orange','darkred']
    for i in range(0,4):
        for j in range(0,4):
            pixel = i*4 + j
            if pixel < 10:
                pixel_label = 'Pixel'+'0'+str(pixel)
            else:
                pixel_label = 'Pixel'+str(pixel)
            for iteration in range(0,len(head_dict[pixel_label][energy])):
                axes[i][j].plot(ThreshV,head_dict[pixel_label][energy][iteration],color = colors[iteration],marker ='.',linewidth=0.5,label='sweep'+str(iteration))
            axes[i][j].plot(ThreshV,np.average(head_dict[pixel_label][energy],axis=0),'ko-')
            axes[i][j].annotate(pixel_label,xy=(1,10),fontsize=12)
            if i == 3 and j == 3:
                axes[i][j].legend(loc='upper right')
            
    big_ax = fig.add_subplot(111)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_color('none')
    big_ax.spines['bottom'].set_color('none')
    big_ax.spines['left'].set_color('none')
    big_ax.spines['right'].set_color('none')
    big_ax.set_ylabel('Counts',fontsize=16)
    big_ax.set_xlabel('Threshold Voltage (V)',fontsize=16)
    big_ax.set_title(title_string,fontsize=16)
    fig.tight_layout()
    fig.savefig(plot_name)
    plt.close(fig)
    
def find_oppt_CDFs(basedir,lower_orbit,upper_orbit,unit):
    
    file_pattern_str  = '%s_%s_engr_l1_sigthrsh.*v1.[23].*.cdf' % serial_num_to_sc(unit)
    file_pattern      = re.compile(file_pattern_str)
    dir_pattern       = re.compile('.*(/\d{4}_201\d).*')
    matching_files    = []
    
    for root, dirs, files in os.walk(basedir):
        for name in files:
            if file_pattern.match(name):
                #get current full path and the parent directory
                curr_path        = root.replace('\\','/')
                #get the orbit_number
                if dir_pattern.match(curr_path):
                    orbit_num_str = get_orb_number(curr_path)
                    #test to see if orbit num is in range       
                    if int(orbit_num_str) >= lower_orbit and int(orbit_num_str) <= upper_orbit:
                        matching_files.append(curr_path+'/'+name)
    
    return matching_files
    
def get_orb_number(file_specification):
    pattern       = re.compile('.*(/\d{4}_201\d).*')
    return pattern.split(file_specification)[1][1:5]