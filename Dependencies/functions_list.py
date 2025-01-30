## Defining some global variables
global omega
root_i = np.sqrt(1j)
omega = [1,root_i,1j,1j*root_i,-1,-root_i,-1j,-1j*root_i]

#A function to calculate n'th power of π/4 (where n is value from 0 to 7)
def expfunc(arr):
    b = np.array([omega[i] for i in arr])
    return b

### Creating the Polynomial array of Circuit
def create_poly(qc, n: int):
    instructions = [(instruction.operation.name,
                    [qc.find_bit(q).index for q in instruction.qubits],
                     instruction.operation.params)
                    for index, instruction in enumerate(qc.data)]

    wire_array = [[i] for i in range(n)]
    # t = n
    max_new_var = n # or the last variable name
    terms = []
    phases = []
    #Each term now also holds a weight, alongside the index of each variable in our polynomial. In form: [weight,[*vars]]

    for entry in instructions:
        gate = entry[0]
        elements = entry[1]
        if gate == 'h':
                wire_array[elements[0]].append(max_new_var)
                #print(wire_array[elements[0]])
                max_new_var += 1
                terms.append([4,[wire_array[elements[0]][-1],wire_array[elements[0]][-2]]]) #Weight for H = 4

        elif gate in ['z','cz','ccz']:
            terms.append([4,[wire_array[j][-1] for j in elements]]) #Weight for Z, CZ, CCZ = 4; only additional thing to note is that CZ and CCZ are multi-qubit gates, so they will give higher degree terms

        elif gate == 's':
            terms.append([2,[wire_array[j][-1] for j in elements]]) #Weight for S = 2

        elif gate == 't':
            terms.append([1,[wire_array[j][-1] for j in elements]]) #Weight for T = 1

        elif gate == 'sdg':
            terms.append([6,[wire_array[j][-1] for j in elements]]) #Weight for Sdg = 6 = (-2) mod 8

        elif gate == 'tdg':
            terms.append([7,[wire_array[j][-1] for j in elements]]) #Weight for Tdg = 7 = (-1) mod 8

        elif gate in ['rz','crz']:
            theta = (entry[2])[0]
            phases.append([theta,[wire_array[j][-1] for j in elements]])

    for term in terms:
        term[-1].sort()
    return terms, wire_array, max_new_var, phases


# function to evaluate polynomial equation
def eval_f(terms,x): # x is a binary array
    val_out = 0
    for term in terms:
        weight = term[0]
        indices = term[1]
        v = bool(1) #The inputs remain boolean. Hence, the products of the variables will remain integers.
        for j in indices:
            v &= x[j]
        val_out = int(val_out + weight*int(v))%8 #Ensuring all operations are done modulo 8, as integer operations.
    return val_out


def truthtable(terms, n, t, initial_state, ivs, np):
    # x_range = number of x values possible, given initial_state
    x_range = 2**(t-n)
    tt = np.empty(x_range,dtype=np.uint8)
    #pt = np.empty(x_range,dtype=np.float128)
    # x = x0 x1 x2 x3 x4 ...
    x = np.empty(t, dtype='bool')
    # y = all variables - input variables == indexes of x2 x3 x4
    y = [i for i in range(t) if i not in ivs]
    # filling input variables value
    tmp = 0
    for j in ivs:
        x[j] = initial_state[tmp]
        tmp += 1
    for i in range(x_range):
        # i = x2x3x4
        # filling other varibles value
        y_bin = bin(i)[2:].zfill(t-n)
        tmp = 0
        for j in y:
            x[j] = bool(int(y_bin[tmp]))
            tmp += 1
        tt[i] = eval_f(terms,x)
    return tt

#Note that we only have to do phase table calculations - to calculate phase factors:
def phasetable(phases, n, t, initial_state, ivs, np):
    x_range = 2**(t-n) # x_range = number of x values possible, given initial_state
    pt = np.ones(x_range,dtype=complex) #This array will store the phase factors (Exp_itheta)

    for term in phases:
        theta = term[0]
        indices = term[1]
        indices.sort()

        if len(indices)==1: #Case for RZ gate

            ind = indices[0]
            if ind < n: #This indicates that the concerned variable is already initialised.

                if initial_state[ind] == 1:

                    pt *= np.exp(1j*theta) #a universal global phase is added, which ultimately does not make a difference but is included for the sake of completeness
            else: #when concerned variable is not initialised and is an intermediate/final variable
                val = np.exp(1j*theta)

                r = t-n
                block_size = 2**(t-ind-1)
                repeat_period = 2**(t-ind)
                for start in range(block_size, x_range, repeat_period):
                    for k in range(block_size):
                        pt[start + k] *= val

        else: #Case for CRZ gate
            ind1 = indices[0]
            ind2 = indices[1]

            if ind1 < n: # if atleast one of the variables is already pre-defined
                if initial_state[ind1] == 1: #this means that we only have to care about ind2 if ind1 index initialised to 1
                    if ind2 < n: # case where ind1 and ind2 are both initialised
                        if initial_state[ind2] == 1:
                            pt *= np.exp(1j*theta)
                    else: #case where ind2 variable is not initialised - equivalent to RZ but on 2nd qubit only.
                        val = np.exp(1j*theta)
                        r = t-n
                        block_size = 2**(t-ind2-1)
                        repeat_period = 2**(t-ind2)
                        for start in range(block_size, x_range, repeat_period):
                            for k in range(block_size):
                                pt[start + k] *= val

          
            else: # if both variables are uninitialised: we go through a nested loop so as to access only states where both i1 and i2 are 1. 
                val = np.exp(1j*theta)
                r = t-n
                block_size_1 = 2**(t-ind1-1)
                repeat_period_1 = 2**(t-ind1)
                block_size_2 = 2**(t-ind2-1)
                repeat_period_2 = 2**(t-ind2)
                for start_1 in range(block_size_1, x_range, repeat_period_1):
                    for start_2 in range(start_1 + block_size_2, repeat_period_1, repeat_period_2):
                        for k in range(block_size_2):
                            pt[start_2+k] *= val
                    
    return pt

def statevector_(ttb, pt, n, t, ovs, np, starting_index=0):
    group_size = 2**(t-n) # == size of ttb
    exp_vals = expfunc(ttb) # == exp(i pi/4 * f(x))
    exp_vals = exp_vals*pt # == Multiply by exp(i theta) for any phase theta to be added

    s = np.zeros(2**n,dtype=complex) 
    for k in range(group_size):
        chosenbits = "".join([(bin(k)[2:].zfill(t))[j] for j in ovs])
        chosen_int = int(chosenbits,2)
        s[chosen_int] += (exp_vals[k])

    stvector = s / (2**.5)**(t-n) #Normalization
    return stvector
