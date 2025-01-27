### Creating the Polynomial array of Circuit
def create_poly(qc, n: int):    
    instructions = [(instruction.operation.name,
                    [qc.find_bit(q).index for q in instruction.qubits]) 
                    for index, instruction in enumerate(qc.data)]

    wire_array = [[i] for i in range(n)]
    # t = n
    max_new_var = n # or the last variable name
    terms = [] 
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


    for term in terms:
        term[-1].sort()
    return terms, wire_array, max_new_var


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


def statevector_(ttb, n, t, ovs, np, starting_index=0):
    group_size = 2**(t-n) # == size of ttb 
    # len(ovs) == n
    s = np.zeros(2**len(ovs),dtype=complex)
    s_ldic = dict()
    for k in range(0,group_size): # Going through each value
        t_val = ttb[k] #Check truth value for each element
        chosenbits = "".join([ ( bin(k)[2:].zfill((t)) )[j] for j in ovs ]) #Choosing the variables which are corresponing to the output. 
        chosen_int = int(chosenbits,2) #Integer value corresponding to chosen variables.  

        #try: s_ldic[chosen_int][t_val]+=1 #If array has been created already, just update it
        #except KeyError: #If the chosen variables have not been chosen before, define a new element corresponding to that combo - and then update the array
        try: s_ldic[chosen_int][t_val]+=1 #If array has been created already, just update it
        except KeyError: #If the chosen variables have not been chosen before, define a new element corresponding to that combo - and then update the array
            s_ldic[chosen_int] = np.array([0,0,0,0,0,0,0,0])
            s_ldic[chosen_int][t_val]+=1

    for k in s_ldic:
        s_ldic[k] = (s_ldic[k][0] - s_ldic[k][4]) + (s_ldic[k][1] - s_ldic[k][5])/np.sqrt(2)- (s_ldic[k][3] - s_ldic[k][7])/np.sqrt(2) + (1j)*((s_ldic[k][2] - s_ldic[k][6]) + (s_ldic[k][1] - s_ldic[k][5])/np.sqrt(2)+ (s_ldic[k][3] - s_ldic[k][7])/np.sqrt(2) ) #Hardcoded the computation of FFT[1] of the array
        s[k] = s_ldic[k] 
        #s_ldic[k] = (s_ldic[k][0] - s_ldic[k][4]) + (s_ldic[k][1] - s_ldic[k][5])/np.sqrt(2)- (s_ldic[k][3] - s_ldic[k][7])/np.sqrt(2) + (1j)*((s_ldic[k][2] - s_ldic[k][6]) + (s_ldic[k][1] - s_ldic[k][5])/np.sqrt(2)+ (s_ldic[k][3] - s_ldic[k][7])/np.sqrt(2) ) #Hardcoded the computation of FFT[1] of the array
        #s[k] = s_ldic[k] 
    stvector = s / (2**.5)**(t-n) #Normalization
    return stvector


