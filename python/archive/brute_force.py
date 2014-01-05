from cantilever_divingboard import *

def brute_force():
    n = 10 # Number of points to check in each dimension

    # ASML
    # 0.45 micron resolution
    # ??? alignment resolution

    f_min = 1e2
    f_max = 100e3
    bandwidth = 100e3

    # Others
    doping = logspace(16, 19, n)
    bias = linspace(0., 10., n)

    # Cantilever dimensions
    length = linspace(10e-6, 100e-6, n)
    width = linspace(0.45e-6, 10e-6, n) # Min = ASML resolution for frontside DRIE
    thickness = linspace(100e-9, 10e-6, n) # Min = stability to survive water entry?

    # For each iteration we'll check that the piezo dimensions are valid first, so match upper limit to length/width/thickness
    piezo_length = linspace(0.45e-6, max(length), n)
    piezo_width = linspace(0.45e-6, max(width), n)
    piezo_gap = linspace(0.45e-6, 10e-6, n)
    piezo_thickness = linspace(100e-9, max(thickness), n)

    # Define limits of iteration
    ii_max = len(doping)
    jj_max = len(bias)
    kk_max = len(length)
    ll_max = len(width)
    mm_max = len(thickness)
    nn_max = len(piezo_length)
    oo_max = len(piezo_width)
    pp_max = len(piezo_gap)
    qq_max = len(piezo_thickness)

    # results = zeros((ii_max, jj_max, kk_max, ll_max, mm_max, nn_max, oo_max, pp_max, qq_max), dtype = Float)

    min_resolution = 10000.

    for ii in range(0, ii_max-1):
        for jj in range(0, jj_max-1):
        
            for kk in range(0, kk_max-1):
                for ll in range(0, ll_max-1):
                    for mm in range(0, mm_max-1):
                
                        for nn in range(0, nn_max-1):
                            for oo in range(0, oo_max-1):
                                for pp in range(0, pp_max-1):
                                    for qq in range(0, qq_max-1):
                                    
                                    
                                        doping_level = doping[ii]
                                        bias_voltage = bias[jj]
                                        l = length[kk]
                                        w = length[ll]
                                        t = length[mm]
                                        l_pr = length[nn]
                                        w_pr = length[oo]
                                        w_gap = length[pp]
                                        t_pr = length[qq]
                                    
                                        # Check dimensions are okay, if not skip onto the next for loop iteration
                                        if (l_pr >= l) or (w_pr >= w) or (t_pr >= t) or (bias_voltage == 0):
                                            continue                                    
                                                    
                                        print bias_voltage                            
                                        c1 = cantilever_divingboard(l, w, t, l_pr, w_pr, w_gap, t_pr, V_bias = bias_voltage, N = doping_level)
                                        
                                        # print 'Cantilever L/W/T: %f %f %f' % (c1.l*1e6, c1.w*1e6, c1.t*1e6)
                                        # print 'Piezo L/W/T: %f %f %f' % (c1.l_pr*1e6, c1.w_pr*1e6, c1.t_pr*1e6)
                                        # print 'Doping: %f' % log10(c1.N)
                                        # print 'Resistance: %f' % c1.resistance()
                                        # print "Integrated noise: %f" % (c1.integrated_noise()*1e8)
                                        # print "Sensitivity: %f" % c1.force_sensitivity()
                                        # print 'Force resolution: %f' % (c1.force_resolution()*1e12)
                                        # print 'Natural freq: %f' % c1.omega_vacuum_hz()
                                        # print '\n'
                                        if (c1.omega_vacuum_hz() > bandwidth):
                                            # results[ii, jj, kk, ll, mm, nn, oo, pp, qq] = c.force_resolution()
                                            # print c1.force_resolution()*1e12
                                            if (c1.force_resolution() < min_resolution):
                                                min_resolution = c1.force_resolution()

    return min_resolution
    
mres = brute_force()
print (mres*1e12)