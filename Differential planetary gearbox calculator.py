'''
This program is intended to generate teeth combinations for
differential planetary gearboxes with sun1 as input, ring1 fixed and 
ring2 as output. For more information visit:
https://juangg-projects.blogspot.ca/2018/02/split-ring-compound-epicyclicplanetary.html
It just tries a bunch of possible combiantions and prints the closest
ones within the parameters and restrictions specified.
Use it at your own risk.

Equations that rule this gearbox:

 gr = 1/((1-(zr1*zp2)/(zr2*zp1))/(1+(zr1/zs1)))
 zr = zs+2*zp
 zp2 = zp1*k
 m1*(zp1+zs1) = m2*(zp2+zs2)
 (zs+zr)%np = 0
 zr2 = zr1*k+np

'''
import pandas as pd
import numpy as np
import math

def max_planets(zs, zp):
    return math.floor(math.pi * (1 + zs / zp))

#Welcome and instructions:
print('''************************************************
*   Differential planetary gearbox calculator  *
*      By Juan Gg                              *
************************************************''')

#Input parameters:
try:
    print('>>> Enter parameters:')
    min_od = float(input(' Min Outside diameter: '))
    max_od = float(input(' Max Outside diameter: '))
    mtr = 40 #float(input(' Min target ratio: '))
    Mtr = 60 #float(input(' Max target ratio: '))
    mt = 11 #int(input(' Min number of teeth: '))
    mm = 1 #float(input(' Min teeth module: '))
    mp = 3 #int(input(' Min number of planets: '))
    Mp = 8 #int(input(' Max number of planets: '))
    if mp > Mp:     #Min number of planets can not be bigger than max number.
        print('  Wrong input!   Try again')


    if input(' Proceed to calculation? (y/n): ').lower() == 'y':
        print('')
        print(' Please wait, it may take several minutes...')
        print('Results: ',end="")
        solutions = []
        r = 0
        # Precompute global max for zr1 based on od and mm
        max_zr1 = math.floor(min_od / mm - 2)
        if max_zr1 < 3 * mt:  # Impossible if max_zr1 too small
            print('No feasible solutions: od too small for mm and mt.')
        else:
            for od in np.arange (min_od,max_od,.1):
                # Precompute global max for zr1 based on od and mm
                max_zr1 = math.floor(od / mm - 2)
                for np in range(mp, Mp + 1):
                    # Per-np bounds for zs1m and zr1m
                    min_zs1m = math.ceil(mt / np)
                    max_zs1m = math.floor((max_zr1 - 2 * mt) / np)  # zr1 = zs1 + 2*zp1 >= zs1 + 2*mt
                    if min_zs1m > max_zs1m:
                        continue  # No feasible zs1 for this np
                    
                    max_zr1m = math.floor(max_zr1 / np)
                    min_delta_m = math.ceil(2 * mt / np)  # Minimum for zr1m - zs1m to make zp1 >= mt
                    
                    k = 0  # Reset k, but we'll limit it per inner config
                    for zs1m in range(min_zs1m, max_zs1m + 1):
                        min_zr1m = max(zs1m + min_delta_m, zs1m + 1)  # zr1 > zs1
                        for zr1m in range(min_zr1m, max_zr1m + 1):
                            diff_m = zr1m - zs1m
                            # Early prune: Ensure zp1 will be integer (np * diff_m must be even)
                            if (np * diff_m) % 2 != 0:
                                continue
                            
                            # Now compute integers directly
                            zs1 = np * zs1m
                            zr1 = np * zr1m
                            zp1 = (zr1 - zs1) // 2  # Integer division, safe due to modulo check
                            
                            # Early module check for stage 1
                            m1 = od / (zr1 + 2)
                            if m1 < mm:
                                continue  # Too many teeth for od
                            
                            # Early planet feasibility check
                            if np > max_planets(zs1, zp1):
                                continue
                            
                            # Now bound k using simplified gr equation
                            # gr = R * (1 + k * (zr1 / np)), where R = 1 + zr1 / zs1
                            R = 1 + zr1 / zs1
                            if R <= 0:  # Degenerate
                                continue
                            
                            # Solve for k range that hits [mtr, Mtr]
                            k_factor = zr1 / np
                            k_min = max(0.0, (mtr / R - 1) / k_factor) if k_factor > 0 else 0.0
                            k_max = (Mtr / R - 1) / k_factor if k_factor > 0 else 0.0
                            if k_min > k_max + 1e-6:  # No overlap
                                continue
                            
                            # Loop k only in feasible range, step 0.1
                            k_start = math.ceil(k_min * 10) / 10.0
                            k_end = math.floor(k_max * 10) / 10.0
                            num_k_steps = int((k_end - k_start) * 10) + 1
                            if num_k_steps <= 0:
                                continue
                            
                            for i in range(num_k_steps + 1):
                                k = k_start + i * 0.1
                                
                                # Compute stage 2 teeth (as floats first)
                                zr2 = k * zr1 + np
                                zp2 = k * zp1
                                zs2 = k * zs1 + np  # Derived: zs2 = zr2 - 2*zp2 = k*zs1 + np
                                
                                # Module for stage 2 (using the relation)
                                denom = k * (zp1 + zs1) + np
                                if denom <= 0:
                                    continue
                                m2 = m1 * (zp1 + zs1) / denom
                                
                                # Verify integers with tolerance (as original)
                                if (zr1 >= mt and zp1 >= mt and zs1 >= mt and
                                    zr2 >= mt and zp2 >= mt and zs2 >= mt and
                                    m1 >= mm and m2 >= mm and
                                    abs(round(zr2) - zr2) < 0.01 and abs(round(zp2) - zp2) < 0.01 and
                                    abs(round(zs2) - zs2) < 0.01):
                                    
                                    # Round to integers
                                    zr1_r = round(zr1)  # Already integer, but for consistency
                                    zp1_r = round(zp1)
                                    zs1_r = round(zs1)
                                    zr2_r = round(zr2)
                                    zp2_r = round(zp2)
                                    zs2_r = round(zs2)
                                    
                                    # Gear ratio
                                    gr = 1 / ((1 - (zr1_r * zp2_r) / (zr2_r * zp1_r)) / (1 + (zr1_r / zs1_r)))
                                    gr = round(gr, 6)
                                    grsp1 = round(zp1_r / zs1_r, 1)
                                    
                                    # Final checks (as original)
                                    if mtr <= gr <= Mtr and np <= max_planets(zs1_r, zp1_r) and np <= max_planets(zs2_r, zp2_r):
                                        r += 1
                                        print(r, end=",")
                                        
                                        # Store (as original, using rounded values)
                                        solutions.append({
                                            'GR': round(gr, 1), 'np': np, 'zr1': zr1_r, 'zp1': zp1_r, 'zs1': zs1_r,
                                            'grsp1': grsp1, 'm1': round(m1, 3), 'zr2': zr2_r, 'zp2': zp2_r,
                                            'zs2': zs2_r, 'm2': round(m2, 3), 'zs2id': round(m2 * (zs2_r - 1.25) - 8, 1)
                                        })
                    
                    # End k outer (but now it's inner; no need for while k<=100)
        
        print('Calculation finished.')
       # Create DataFrame from the list
        solutions_df = pd.DataFrame(solutions)
        # choose the subset based on the sun gear diameter 
        solutions_df = solutions_df[(solutions_df['zs2id'] >= 20)]
        # eliminate scale duplicates
        solutions_df.sort_values(by='m1')
        # Step 2: Define subset columns (all except 'm1')
        subset_cols = [col for col in solutions_df.columns if col != 'm1' and col != 'm2' and col != 'zs2id']
        # Step 3: Drop duplicates on subset, keep first (min m1)
        solutions_df = solutions_df.drop_duplicates(subset=subset_cols).reset_index(drop=True)
        
        # sort the values by the gear ratio of the input on gear to the planets
        solutions_df = solutions_df.sort_values(by='grsp1')
        print(solutions_df.to_string())
except Exception as e:
    import traceback
    print(f"Error type: {type(e).__name__}")
    print(f"Error message: {str(e)}")
    traceback.print_exc()  # This prints the full stack trace with line numbers