from affine import *
from common import *
from cycle import *

def is_commute(f, g):
    a_1 = f.a
    b_1 = f.b
    a_2 = g.a
    b_2 = g.b
    if ((a_1-1)*b_2 - (a_2-1)*b_1)%P == 0:
        return True
    else:
        return False
    
def check_req_1(f_list, g_list):
    for i in range(l_h):
        for j in range(-(J-1), J):
            f_idx = i
            g_idx = (-i + j) % l_h
            if not is_commute(f_list[f_idx], g_list[g_idx]):
                return False
    return True
            
def check_inner_commute(f_list):
    flag = False
    for i in range(l_h):
        for j in range(i+1, l_h):
            if  is_commute(f_list[i], f_list[j]) == False:
                flag = True
    return flag
def check_condition_a(f_list, g_list):
    return check_req_1(f_list, g_list)

def check_condition_b(f_list, g_list):
    return check_inner_commute(f_list) and check_inner_commute(g_list)

def check_condition_c(f_list, g_list, cols_list):
    h_x = make_h_x_functions(f_list, g_list)
    h_z = make_h_z_functions(f_list, g_list)
    cbcs = []
    for cols in cols_list:
        if not is_utcbc(cols):
            cycle_function_x = [AffineFunc(1, 0), AffineFunc(1, 0)]
            cycle_function_z = [AffineFunc(1, 0), AffineFunc(1, 0)]
            cycles = [make_cycle(cols, 0), make_cycle(cols, 1)]
            for i in range(int(len(cycles[0])/2)):
                cycle_function_x[0] = composite_affine(cycle_function_x[0], h_x[cycles[0][2*i][0]][cycles[0][2*i][0]])
                cycle_function_x[0] = composite_affine(cycle_function_x[0], h_x[cycles[0][2*i+1][0]][cycles[0][2*i+1][0]].inv())
                cycle_function_x[1] = composite_affine(cycle_function_x[1], h_x[cycles[1][2*i][0]][cycles[1][2*i][0]])
                cycle_function_x[1] = composite_affine(cycle_function_x[1], h_x[cycles[1][2*i+1][0]][cycles[1][2*i+1][0]].inv())
                cycle_function_z[0] = composite_affine(cycle_function_z[0], h_z[cycles[0][2*i][0]][cycles[0][2*i][0]])
                cycle_function_z[0] = composite_affine(cycle_function_z[0], h_z[cycles[0][2*i+1][0]][cycles[0][2*i+1][0]].inv())
                cycle_function_z[1] = composite_affine(cycle_function_z[1], h_z[cycles[1][2*i][0]][cycles[1][2*i][0]])
                cycle_function_z[1] = composite_affine(cycle_function_z[1], h_z[cycles[1][2*i+1][0]][cycles[1][2*i+1][0]].inv())
            if is_cbc(cycle_function_x[0]) or is_cbc(cycle_function_x[1]) or is_cbc(cycle_function_z[0]) or is_cbc(cycle_function_z[1]):
                cbcs.append(cols)
    if cbcs:
        return cbcs
    else:
        return True

def is_cbc(affine_func):
    if affine_func.a == 1:
        if affine_func.b == 0:
            return True
        else:
            return False
    else:
        if math.gcd(affine_func.a-1, P):
            return True
        else:
            return False

def check_conditions(f_list, g_list, cols_list):
    if check_condition_a(f_list, g_list) and check_condition_b(f_list, g_list) and check_condition_c(f_list, g_list, cols_list):
        return True
    return False