### Creating the Polynomial array of Circuit
def create_poly(qc, n: int):    
    instructions = [(instruction.operation.name,
                    [qc.find_bit(q).index for q in instruction.qubits]) 
                    for index, instruction in enumerate(qc.data)]

    wire_array = [[i] for i in range(n)]
    t = n
    terms = []

    for entry in instructions:
        gate = entry[0]
        elements = entry[1]

    if gate == 'h':
            wire_array[elements[0]].append(max_new_var)
            #print(wire_array[elements[0]])
            max_new_var+=1
            terms.append([4,[wire_array[elements[0]][-1],wire_array[elements[0]][-2]]])


        elif gate in ['z','cz','ccz']:
            terms.append([4,[wire_array[j][-1] for j in elements]])

        elif gate == 's':
            terms.append([2,[wire_array[j][-1] for j in elements]])

        elif gate == 't':
            terms.append([1,[wire_array[j][-1] for j in elements]])
        
        elif gate == 'sdg':
            terms.append([6,[wire_array[j][-1] for j in elements]])

        elif gate == 'tdg':
            terms.append([7,[wire_array[j][-1] for j in elements]])


    for term in terms:
        term[-1].sort()
    return terms, wire_array


# function to evaluate polynomial equation
def eval_f(terms,x): # x is a binary array 
    val_out = bool(0)
    for term in terms:
        v = bool(1)
        for j in term:
            v &= x[j]
        val_out += (weight*int(v))%8
        val_out %=8
    return val_out

def truthtable(terms, n, t, initial_state, ivs, np):
    # x_range = number of x values possible, given initial_state
    x_range = 2**(t-n) 
    tt = np.empty(x_range,dtype='bool')
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


def statevector_(ttb, n, t, ovs, np):
    group_size = 2**(t-n) # == size of ttb 
    # len(ovs) == n
    s = np.zeros(2**len(ovs),dtype=complex)
    s_ldic = dict()
   for k in range(starting_index,starting_index+group_size): # Going through each value
        t_val = ttb[k,0]
        chosenbits = "".join([(bin(k)[2:].zfill((max_new_var)))[j] for j in ovs])
        chosen_int = int(chosenbits,2)

        try: s_ldic[chosen_int][t_val]+=1
        except KeyError:    
            s_ldic[chosen_int] = np.array([0,0,0,0,0,0,0,0])
            s_ldic[chosen_int][t_val]+=1

    for k in s_ldic:
      s[k] = (s_ldic[k][0] - s_ldic[k][4]) + (s_ldic[k][1] - s_ldic[k][5])/np.sqrt(2)- (s_ldic[k][3] - s_ldic[k][7])/np.sqrt(2) + (1j)*((s_ldic[k][2] - s_ldic[k][6]) + (s_ldic[k][1] - s_ldic[k][5])/np.sqrt(2)+ (s_ldic[k][3] - s_ldic[k][7])/np.sqrt(2) )
    stvector = s / (2**.5)**(t-n)
    return stvector


