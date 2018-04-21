#!/usr/bin/env python
import fileinput
import re


#Dirty search-replace loop without regex. But it works!
#Would probably benefit from some clean regex code > to do list
 
for line in fileinput.input():
    line = re.sub('QUERY___rdgB_hit_','', line.rstrip())
    line = re.sub('QUERY___PLC_hit_','', line.rstrip())
    line = re.sub('QUERY___Gprk2_hit_','', line.rstrip())
    line = re.sub('QUERY___r_opsin_hit_','', line.rstrip())
    line = re.sub('QUERY___trp_hit_','', line.rstrip())
    line = re.sub('QUERY___Gq_beta_hit_','', line.rstrip())
    line = re.sub('QUERY___Gq_gamma_hit_','', line.rstrip())
    line = re.sub('QUERY___Gq_alpha_hit_','', line.rstrip())
    line = re.sub('QUERY___Gprk1_hit_','', line.rstrip())
    line = re.sub('QUERY___DAGK_hit_','', line.rstrip())
    line = re.sub('QUERY___Arr_hit_','', line.rstrip())
    line = re.sub('QUERY___rdgC_hit_','', line.rstrip())
    line = re.sub('QUERY___PKC_hit_','', line.rstrip())
    line = re.sub('_ORF1','', line.rstrip())
    line = re.sub('_ORF2','', line.rstrip())
    line = re.sub('_ORF3','', line.rstrip())
    line = re.sub('_ORF4','', line.rstrip())
    line = re.sub('_ORF5','', line.rstrip())
    line = re.sub('_ORF6','', line.rstrip())
    line = re.sub('_ORF7','', line.rstrip())
    line = re.sub('_ORF8','', line.rstrip())
    line = re.sub('_ORF9','', line.rstrip())
    line = re.sub('_ORF10','', line.rstrip())
    line = re.sub('_ORF11','', line.rstrip())
    line = re.sub('_ORF12','', line.rstrip())
    line = re.sub('\\|','_', line.rstrip())

    print(line)
