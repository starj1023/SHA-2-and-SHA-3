from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, Swap, Z, T, Tdag, S, Tdagger, Sdag
from projectq.backends import CircuitDrawer, ResourceCounter, CommandPrinter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control, Dagger
from math import floor, ceil, log2


def quantum_and(a, b, c, ancilla):

    H | c
    CNOT | (b, ancilla)
    CNOT | (c, a)
    CNOT | (c, b)
    CNOT | (a, ancilla)
    Tdag | a
    Tdag | b
    T | c
    T | ancilla
    CNOT | (a, ancilla)
    CNOT | (c, b)
    CNOT | (c, a)
    CNOT | (b, ancilla)
    H | c
    S | c

def quantum_and_dag(a, b, c):

    H | b
    H | c
    Measure | c
    if(int(c) or resource_check == 1):
        CNOT | (a,b)
        X | c
    H | b

def toffoli_gate(a, b, c, ancilla = 0, mode=True):

    if(TD == 1):
        Tdag | a
        Tdag | b
        H | c
        CNOT | (c, a)
        T | a
        CNOT | (b, c)
        CNOT | (b, a)
        T | c
        Tdag | a
        CNOT | (b, c)
        CNOT | (c, a)
        T | a
        Tdag | c
        CNOT | (b, a)
        H | c
    elif (TD == 2):
        if mode:
            quantum_and(a,b,c,ancilla)
        else:
            quantum_and_dag(a,b,c)
    else:
        Toffoli | (a, b, c)


def w(n): # for draper
    return n - sum(int(floor(n / (pow(2, i)))) for i in range(1, int(log2(n)) + 1))

def l(n, t): # for draper
    return int(floor(n / (pow(2, t))))

def outDraper(eng,a,b, ancillas, output = 0):
    n = len(a)

    if output != 0:
        z = output[:n]
    else:
        z = eng.allocate_qureg(n)
    length = n-1-w(n-1)-floor(log2(n-1))
    ancilla = ancillas[:length]
    and_len = n-1
    and_idx = 0
    tmp_len = (n-1)-w(n-1)
    if TD == 2:
        and_ancilla = ancillas[length:length+and_len]
        tmp_ancilla = ancillas[length+and_len:]
    else:
        and_ancilla = 0
        tmp_ancilla = 0
    tmp_idx = 0

    # Init round
    for i in range(n-1):
        if TD == 2:
            toffoli_gate(a[i], b[i], z[i + 1], and_ancilla[and_idx]) # and
        else:
            toffoli_gate(a[i], b[i], z[i + 1])
        and_idx += 1
    for i in range(1,n-1):
        CNOT | (a[i], b[i])
    CNOT | (a[n-1], z[n-1])
    CNOT | (b[n-1], z[n-1])

    # P-round
    #print("P-rounds")
    idx = 0 # ancilla idx
    tmp = 0
    and_idx = 0
    for t in range(1, int(log2(n-1))):
        pre = tmp
        for m in range(1, l(n-1, t)):
            if t == 1:
                if TD == 2:
                    toffoli_gate(b[2 * m], b[2 * m + 1], ancilla[idx], and_ancilla[and_idx])
                else:
                    toffoli_gate(b[2 * m], b[2 * m + 1], ancilla[idx])
            else:
                if TD == 2:
                    toffoli_gate(ancilla[pre - 1 + 2 * m], ancilla[pre - 1 + 2 * m + 1], ancilla[idx], and_ancilla[and_idx])
                else:
                    toffoli_gate(ancilla[pre - 1 + 2 * m], ancilla[pre - 1 + 2 * m + 1], ancilla[idx])
            if m == 1:
                tmp = idx
            idx += 1
            and_idx += 1

    # G-round
    #print("G-rounds")
    pre = 0  # The number of cumulative p(t-1)
    idx = 0  # ancilla idx
    for t in range(1, int(log2(n-1)) + 1):
        for m in range(l(n-1, t)):
            if t == 1:
                if TD == 2:
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], b[2 * m + 1], tmp_ancilla[tmp_idx], and_ancilla[and_idx])
                    CNOT | (tmp_ancilla[tmp_idx], z[int(pow(2, t) * (m + 1))])
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], b[2 * m + 1], tmp_ancilla[tmp_idx], and_ancilla[and_idx], False)
                else:
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], b[2 * m + 1], z[int(pow(2, t) * (m + 1))])

            else:
                if TD == 2:
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], ancilla[idx+2*m], tmp_ancilla[tmp_idx], and_ancilla[and_idx])
                    CNOT | (tmp_ancilla[tmp_idx], z[int(pow(2, t) * (m + 1))])
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], ancilla[idx+2*m], tmp_ancilla[tmp_idx], and_ancilla[and_idx], False)
                else:
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], ancilla[idx+2*m], z[int(pow(2, t) * (m + 1))])
            and_idx = (and_idx+1) % and_len
            tmp_idx += 1
        if t != 1:
            pre = pre + l(n-1, t-1) -1
            idx = pre

    # C-round
    #print("C-rounds")
    if int(log2(n-1)) - 1 == int(log2(2 * (n-1) / 3)):
        iter = l(n-1, int(log2(n-1)) - 1) - 1
    else:
        iter = 0
    pre = 0
    tmp_idx = 0
    for t in range(int(log2(2 * (n-1) / 3)), 0, -1):
        for m in range(1, l((n-1 - pow(2, t-1)),t)+1):
            if t == 1:
                if TD == 2:
                    toffoli_gate(z[int(pow(2, t) * m)], b[2 * m], tmp_ancilla[tmp_idx], and_ancilla[and_idx])
                    CNOT | (tmp_ancilla[tmp_idx], z[int(pow(2, t) * m + pow(2, t - 1))])
                    toffoli_gate(z[int(pow(2, t) * m)], b[2 * m], tmp_ancilla[tmp_idx], and_ancilla[and_idx], False)
                else:
                    toffoli_gate(z[int(pow(2, t) * m)], b[2 * m], z[int(pow(2, t) * m + pow(2, t - 1))])
            else:
                if m==1:
                    iter += l(n-1, t - 1) - 1
                    pre = length - 1 - iter
                if TD == 2:
                    toffoli_gate(z[int(pow(2, t) * m)], ancilla[pre + 2 * m],tmp_ancilla[tmp_idx], and_ancilla[and_idx])
                    CNOT | (tmp_ancilla[tmp_idx],z[int(pow(2, t) * m + pow(2, t - 1))])
                    toffoli_gate(z[int(pow(2, t) * m)], ancilla[pre + 2 * m],tmp_ancilla[tmp_idx], and_ancilla[and_idx], False)
                else:
                    toffoli_gate(z[int(pow(2, t) * m)], ancilla[pre + 2 * m], z[int(pow(2, t) * m + pow(2, t-1))])

            and_idx = (and_idx + 1) % and_len
            tmp_idx += 1

    # P-inverse round
    #print("P-inv-rounds")
    pre = 0
    iter = l(n-1, int(log2(n-1)) - 1) - 1
    iter2 = 0 # for idx
    idx = 0
    and_idx = 0
    for t in reversed(range(1, int(log2(n-1)))):
        for m in range(1, l(n-1, t)):
            if t == 1:
                if TD == 2:
                    toffoli_gate(b[2 * m], b[2 * m + 1], ancilla[m - t], and_ancilla[and_idx], False)
                else:
                    toffoli_gate(b[2 * m], b[2 * m + 1], ancilla[m - t])
            else:
                if m == 1:
                    iter += l(n-1, t - 1) - 1  # p(t-1) last idx
                    pre = length - iter
                    iter2 += (l(n-1, t) - 1)
                    idx = length - iter2
                if TD == 2:
                    toffoli_gate(ancilla[pre - 1 + 2 * m], ancilla[pre - 1 + 2 * m + 1], ancilla[idx-1+m], and_ancilla[and_idx], False)
                else:
                    toffoli_gate(ancilla[pre - 1 + 2 * m], ancilla[pre - 1 + 2 * m + 1], ancilla[idx-1+m])
            and_idx += 1

    # Last round
    for i in range(n-1):
        CNOT | (b[i], z[i])
    CNOT | (a[0], z[0])
    for i in range(1, n-1):
        CNOT | (a[i], b[i])

    return z

def outDraper_dag(eng, a,b, ancillas, output):
    n = len(a)
    z = output[:n]
    length = n-1-w(n-1)-floor(log2(n-1))
    #ancilla = eng.allocate_qureg(length)
    ancilla = ancillas[:length]
    and_len = n-1
    #and_ancilla = eng.allocate_qureg(and_len)
    and_ancilla = ancillas[length:length+and_len]
    and_idx = 0
    tmp_len = (n-1)-w(n-1)
    #tmp_ancilla = eng.allocate_qureg(tmp_len)
    tmp_ancilla = ancillas[length+and_len:]
    tmp_idx = 0

    # Last round
    for i in range(1, n - 1):
        CNOT | (a[i], b[i])
    CNOT | (a[0], z[0])
    for i in range(n - 1):
        CNOT | (b[i], z[i])

    # P-round
    idx = 0  # ancilla idx
    tmp = 0
    for t in range(1, int(log2(n - 1))):
        pre = tmp
        for m in range(1, l(n - 1, t)):
            if t == 1:
                if TD == 2:
                    toffoli_gate(b[2 * m], b[2 * m + 1], ancilla[idx], and_ancilla[and_idx])  # and
                else:
                    toffoli_gate(b[2 * m], b[2 * m + 1], ancilla[idx])
            else:
                # print(pre - 1 + 2 * m,pre - 1 + 2 * m + 1,idx)
                if TD == 2:
                    toffoli_gate(ancilla[pre - 1 + 2 * m], ancilla[pre - 1 + 2 * m + 1], ancilla[idx],
                             and_ancilla[and_idx])
                else:
                    toffoli_gate(ancilla[pre - 1 + 2 * m], ancilla[pre - 1 + 2 * m + 1], ancilla[idx])
            if m == 1:
                tmp = idx
            idx += 1
            and_idx += 1

    # C-round - reverse
    pre = 0
    tmp_idx = 0
    for t in reversed(range(int(log2(2 * (n - 1) / 3)), 0, -1)):
        idx = pre
        for m in range(1, l((n - 1 - pow(2, t - 1)), t) + 1):
            if t == 1:
                if TD == 2:
                    toffoli_gate(z[int(pow(2, t) * m)], b[2 * m], tmp_ancilla[tmp_idx], and_ancilla[and_idx])
                    CNOT | (tmp_ancilla[tmp_idx], z[int(pow(2, t) * m + pow(2, t - 1))])
                    toffoli_gate(z[int(pow(2, t) * m)], b[2 * m], tmp_ancilla[tmp_idx], and_ancilla[and_idx], False)
                else:
                    toffoli_gate(z[int(pow(2, t) * m)], b[2 * m], z[int(pow(2, t) * m + pow(2, t - 1))])
            else:
                if TD == 2:
                    toffoli_gate(z[int(pow(2, t) * m)], ancilla[idx - 1 + 2 * m], tmp_ancilla[tmp_idx],
                                 and_ancilla[and_idx])
                    CNOT | (tmp_ancilla[tmp_idx], z[int(pow(2, t) * m + pow(2, t - 1))])
                    toffoli_gate(z[int(pow(2, t) * m)], ancilla[idx - 1 + 2 * m], tmp_ancilla[tmp_idx],
                                 and_ancilla[and_idx], False)
                else:
                    toffoli_gate(z[int(pow(2, t) * m)], ancilla[idx - 1 + 2 * m], z[int(pow(2, t) * m + pow(2, t - 1))])
                if m == 1:
                    pre += l(n - 1, t - 1) - 1
            and_idx = (and_idx + 1) % and_len
            tmp_idx += 1

    # G-round - reverse
    pre = 0  # The number of cumulative p(t-1)
    idx_t = int(log2(n - 1))
    iter = 0
    tmp_idx = 0
    for t in reversed(range(1, int(log2(n - 1)) + 1)):
        for m in range(l(n - 1, t)):
            if t == 1:
                if TD == 2:
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], b[2 * m + 1], tmp_ancilla[tmp_idx],
                                 and_ancilla[and_idx])
                    CNOT | (tmp_ancilla[tmp_idx], z[int(pow(2, t) * (m + 1))])
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], b[2 * m + 1], tmp_ancilla[tmp_idx],
                                 and_ancilla[and_idx], False)
                else:
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], b[2 * m + 1], z[int(pow(2, t) * (m + 1))])
            else:
                if m == 0:
                    iter += l(n - 1, idx_t - 1) - 1  # p(t-1) last idx
                    pre = length - iter
                    idx_t -= 1

                if TD == 2:
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], ancilla[pre - 1 + 2 * m + 1], tmp_ancilla[tmp_idx],
                                 and_ancilla[and_idx])
                    CNOT | (tmp_ancilla[tmp_idx], z[int(pow(2, t) * (m + 1))])
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], ancilla[pre - 1 + 2 * m + 1], tmp_ancilla[tmp_idx],
                                 and_ancilla[and_idx], False)
                else:
                    toffoli_gate(z[int(pow(2, t) * m + pow(2, t - 1))], ancilla[pre - 1 + 2 * m + 1], z[int(pow(2, t) * (m + 1))])

            and_idx = (and_idx + 1) % and_len
            tmp_idx += 1

    # P-inverse round
    pre = 0
    iter = l(n - 1, int(log2(n - 1)) - 1) - 1
    iter2 = 0  # for idx
    idx = 0
    and_idx = 0
    for t in reversed(range(1, int(log2(n - 1)))):
        for m in range(1, l(n - 1, t)):
            if t == 1:
                if TD == 2:
                    toffoli_gate(b[2 * m], b[2 * m + 1], ancilla[m - t], and_ancilla[and_idx], False)
                else:
                    toffoli_gate(b[2 * m], b[2 * m + 1], ancilla[m - t])

            else:
                if m == 1:
                    iter += l(n - 1, t - 1) - 1  # p(t-1) last idx
                    pre = length - iter
                    iter2 += (l(n - 1, t) - 1)
                    idx = length - iter2

                if TD == 2:
                    toffoli_gate(ancilla[pre - 1 + 2 * m], ancilla[pre - 1 + 2 * m + 1], ancilla[idx - 1 + m],
                                 and_ancilla[and_idx], False)
                else:
                    toffoli_gate(ancilla[pre - 1 + 2 * m], ancilla[pre - 1 + 2 * m + 1], ancilla[idx - 1 + m])
            and_idx += 1

    and_idx = 0
    # Init round
    CNOT | (b[n - 1], z[n - 1])
    CNOT | (a[n - 1], z[n - 1])
    for i in range(1, n - 1):
        CNOT | (a[i], b[i])
    for i in range(n - 1):
        if TD == 2:
            toffoli_gate(a[i], b[i], z[i + 1], 0, False)
        else:
            toffoli_gate(a[i], b[i], z[i + 1])
        and_idx += 1

    return z

def modQFA(a,b,c):
    CNOT | (a,b)
    CNOT | (b,c)
    CNOT | (a,b)

def modQFA_R(a, b, c):
    CNOT | (a, b)
    CNOT | (b, c)
    CNOT | (a, b)

def QFA(a, b, c, d, ancilla=0):
    CNOT | (a, b)
    CNOT | (a, c)
    if TD == 2:
        toffoli_gate(b, c, d, ancilla)
    else:
        toffoli_gate(b, c, d)
    CNOT | (a, b)
    CNOT | (a, d)
    CNOT | (b, c)

def QFA_R(a, b, c, d):
    CNOT | (b, c)
    CNOT | (a, d)
    CNOT | (a, b)
    if TD == 2:
        toffoli_gate(b, c, d, 0, False)
    else:
        toffoli_gate(b, c, d)
    CNOT | (a, c)
    CNOT | (a, b)

def QHA(a, b, c, ancilla=0):
    if TD == 2:
        toffoli_gate(a, b, c, ancilla)
    else:
        toffoli_gate(a, b, c)
    CNOT | (a, b)

def QHA_R(a, b, c):
    CNOT | (a, b)
    if TD == 2:
        toffoli_gate(a, b, c, 0, False)
    else:
        toffoli_gate(a, b, c)

def initialize_circuit_and_qubits(input):
    n = 32
    w = []
    for i in range(n):
        tmp = []
        for qb in input:
            if len(qb) == n-1:
                if i != 0:
                    tmp.append(qb[i-1])
            else:
                tmp.append(qb[i])
        w.append(tmp)

    return w

global an
an = []
global share1
global share2
share1 = []
share2 = []

def apply_adders_and_handle_carry00(eng, w, input2, ancilla, and_ancilla):
    num = 0
    and_idx = 0
    t = 1
    copy_ancilla = ancilla[-32:]
    idx = 0
    stack = []
    share1 = []
    share2 = []
    flag1 = 1
    flag2 = 1
    while (t):
        n = len(w)
        t = 0
        w_tmp = []

        if (len(w[n - 1]) > 2):
            while (len(w[n - 1]) > 2):
                a = w[n - 1].pop(0)
                b = w[n - 1].pop(0)
                c = w[n - 1].pop(0)
                stack.append([a, b, c, 0, 0, 0])
                modQFA(a, b, c)
                w_tmp.append(c)

            if flag1 == 1:
                flag1 = 0
                stack.append([w_tmp[0],copy_ancilla[idx]])
                CNOT | (w_tmp[0], copy_ancilla[idx])
                share1.append(copy_ancilla[idx])
                idx+=1
                w[n - 1].insert(0, w_tmp[0])
                w[n - 1].insert(-1, w_tmp[1])

            else:
                w[n - 1].extend(w_tmp)

        for i in reversed(range(0, n - 1)):
            w_tmp = []
            w_tmp2 = []
            while (len(w[i]) > 2):
                a = w[i].pop(0)
                b = w[i].pop(0)
                c = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, c, carry])
                if TD == 2:
                    QFA(a, b, c, carry, and_ancilla[and_idx])
                else:
                    QFA(a, b, c, carry)
                and_idx += 1
                w_tmp.append(c)
                w_tmp2.append(carry)

            if flag2 == 1:
                stack.append([w_tmp[0], copy_ancilla[idx]])
                CNOT | (w_tmp[0], copy_ancilla[idx])
                share1.append(copy_ancilla[idx])
                idx += 1
                share2.append(w_tmp2[0])
                w[i].insert(0, w_tmp[0])
                w[i].insert(-1, w_tmp[1])
                w[i + 1].insert(1, w_tmp2[0])
                w[i + 1].insert(-1, w_tmp2[1])

            else:
                w[i].extend(w_tmp)
                w[i + 1].extend(w_tmp2)

        and_idx = 0
        if flag2 == 1:
            flag2 = 0
            share1.reverse()
            share2.reverse()
            input2.append(share1)
            input2.append(share2)

    return w, stack, input2, num

def apply_adders_and_handle_carry0(eng, w, input2, ancilla, and_ancilla=0):
    num = 0
    and_idx = 0
    t = 1
    stack = []
    share1 = []
    share2 = []
    flag1 = 1
    flag2 = 1
    while (t):
        n = len(w)
        t = 0
        w_tmp = []

        if (len(w[n - 1]) > 2):
            while (len(w[n - 1]) > 2):
                a = w[n - 1].pop(0)
                b = w[n - 1].pop(0)
                c = w[n - 1].pop(0)
                stack.append([a, b, c, 0, 0, 0])
                modQFA(a, b, c)
                w_tmp.append(c)

            if flag1 == 1:
                flag1 = 0
                share1.append(w_tmp[0])
                w[n - 1].insert(0, w_tmp[0])
                w[n - 1].insert(-1, w_tmp[1])

            else:
                w[n - 1].extend(w_tmp)

        for i in reversed(range(0, n - 1)):
            w_tmp = []
            w_tmp2 = []
            while (len(w[i]) > 2):
                a = w[i].pop(0)
                b = w[i].pop(0)
                c = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, c, carry])
                if TD == 2:
                    QFA(a, b, c, carry, and_ancilla[and_idx])
                else:
                    QFA(a, b, c, carry)
                and_idx += 1
                w_tmp.append(c)
                w_tmp2.append(carry)

            if flag2 == 1:
                share1.append(w_tmp[0])
                share2.append(w_tmp2[0])
                w[i].insert(0, w_tmp[0])
                w[i].insert(-1, w_tmp[1])
                w[i + 1].insert(1, w_tmp2[0])
                w[i + 1].insert(-1, w_tmp2[1])

            else:
                w[i].extend(w_tmp)
                w[i + 1].extend(w_tmp2)

        and_idx = 0
        if flag2 == 1:
            flag2 = 0
            share1.reverse()
            share2.reverse()
            input2.append(share1)
            input2.append(share2)

    return w, stack, input2, num

def apply_adders_and_handle_carry01(eng, w, stack, ancilla, num, and_ancilla):
    and_idx = 0
    t = 1
    flag = 0
    while (t):
        n = len(w)
        t = 0
        w_tmp = []
        if (len(w[n - 1]) > 2):
            while (len(w[n - 1]) > 2):
                a = w[n - 1].pop(0)
                b = w[n - 1].pop(0)
                c = w[n - 1].pop(0)
                stack.append([a, b, c, 0, 0, 0])
                modQFA(a, b, c)
                w_tmp.append(c)

            w[n - 1].extend(w_tmp)

        for i in reversed(range(0, n - 1)):
            w_tmp = []
            w_tmp2 = []
            while (len(w[i]) > 2):
                a = w[i].pop(0)
                b = w[i].pop(0)
                c = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, c, carry])
                if TD == 2:
                    QFA(a, b, c, carry, and_ancilla[and_idx])
                else:
                    QFA(a, b, c, carry)
                and_idx += 1
                w_tmp.append(c)
                w_tmp2.append(carry)

            if i == 0:
                flag = 1

            if flag == 1 and len(w[i]) == 2:
                a = w[i].pop(0)
                b = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, carry])
                if TD == 2:
                    QHA(a, b, carry, and_ancilla[and_idx])
                else:
                    QHA(a, b, carry)
                and_idx += 1
                w_tmp.append(b)
                w_tmp2.append(carry)

            w[i].extend(w_tmp)
            w[i + 1].extend(w_tmp2)

        and_idx = 0
        for i in range(len(w)):
            if (len(w[i]) > 3):
                t += 1
                flag = 1

    and_idx = 0
    t = 1
    while (t):
        n = len(w)
        t = 0
        w_tmp = []
        if (len(w[n - 1]) > 2):
            while (len(w[n - 1]) > 2):
                a = w[n - 1].pop(0)
                b = w[n - 1].pop(0)
                c = w[n - 1].pop(0)
                stack.append([a, b, c, 0, 0, 0])
                modQFA(a, b, c)
                w_tmp.append(c)

            w[n - 1].extend(w_tmp)

        for i in reversed(range(0, n - 1)):
            w_tmp = []
            w_tmp2 = []
            while (len(w[i]) > 2):
                a = w[i].pop(0)
                b = w[i].pop(0)
                c = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, c, carry])
                if TD == 2:
                    QFA(a, b, c, carry, and_ancilla[and_idx])
                else:
                    QFA(a, b, c, carry)
                and_idx += 1
                w_tmp.append(c)
                w_tmp2.append(carry)

            if len(w[i]) == 2:
                a = w[i].pop(0)
                b = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, carry])
                if TD == 2:
                    QHA(a, b, carry, and_ancilla[and_idx])
                else:
                    QHA(a, b, carry)
                and_idx += 1
                w_tmp.append(b)
                w_tmp2.append(carry)

            w[i].extend(w_tmp)
            w[i + 1].extend(w_tmp2)

        for i in range(len(w)):
            if (len(w[i]) > 2):
                t += 1
    return w, stack

def apply_adders_and_handle_carry1(eng, w, stack, ancilla, num, and_ancilla=0):
    and_idx = 0
    t = 1
    flag = 0
    while (t):
        n = len(w)
        t = 0
        w_tmp = []
        if (len(w[n - 1]) > 2):
            while (len(w[n - 1]) > 2):
                a = w[n - 1].pop(0)
                b = w[n - 1].pop(0)
                c = w[n - 1].pop(0)
                stack.append([a, b, c, 0, 0, 0])
                modQFA(a, b, c)
                w_tmp.append(c)

            w[n - 1].extend(w_tmp)

        for i in reversed(range(0, n - 1)):
            w_tmp = []
            w_tmp2 = []
            while (len(w[i]) > 2):
                a = w[i].pop(0)
                b = w[i].pop(0)
                c = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, c, carry])
                if TD == 2:
                    QFA(a, b, c, carry, and_ancilla[and_idx])
                else:
                    QFA(a, b, c, carry)
                and_idx += 1
                w_tmp.append(c)
                w_tmp2.append(carry)

            if flag == 1 and len(w[i]) == 2:
                a = w[i].pop(0)
                b = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, carry])
                if TD == 2:
                    QHA(a, b, carry, and_ancilla[and_idx])
                else:
                    QHA(a, b, carry)
                and_idx += 1
                w_tmp.append(b)
                w_tmp2.append(carry)

            w[i].extend(w_tmp)
            w[i + 1].extend(w_tmp2)

        and_idx = 0
        for i in range(len(w)):
            if (len(w[i]) > 3):
                t += 1
                flag = 1

    and_idx = 0
    t = 1
    while (t):
        n = len(w)
        t = 0
        w_tmp = []
        if (len(w[n - 1]) > 2):
            while (len(w[n - 1]) > 2):
                a = w[n - 1].pop(0)
                b = w[n - 1].pop(0)
                c = w[n - 1].pop(0)
                stack.append([a, b, c, 0, 0, 0])
                modQFA(a, b, c)
                w_tmp.append(c)

            w[n - 1].extend(w_tmp)

        for i in reversed(range(0, n - 1)):
            w_tmp = []
            w_tmp2 = []
            while (len(w[i]) > 2):
                a = w[i].pop(0)
                b = w[i].pop(0)
                c = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, c, carry])
                if TD == 2:
                    QFA(a, b, c, carry, and_ancilla[and_idx])
                else:
                    QFA(a, b, c, carry)
                and_idx += 1
                w_tmp.append(c)
                w_tmp2.append(carry)

            if len(w[i]) == 2:
                a = w[i].pop(0)
                b = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, carry])
                if TD == 2:
                    QHA(a, b, carry, and_ancilla[and_idx])
                else:
                    QHA(a, b, carry)
                and_idx += 1
                w_tmp.append(b)
                w_tmp2.append(carry)

            w[i].extend(w_tmp)
            w[i + 1].extend(w_tmp2)

        and_idx = 0
        for i in range(len(w)):
            if (len(w[i]) > 2):
                t += 1
    return w, stack

def apply_adders_and_handle_carry(eng, w, ancilla, and_ancilla):
    and_idx = 0
    num = 0
    t = 1
    stack = []
    flag = 0
    while (t):
        n = len(w)
        t = 0
        w_tmp = []
        if (len(w[n - 1]) > 2):
            while (len(w[n - 1]) > 2):
                a = w[n - 1].pop(0)
                b = w[n - 1].pop(0)
                c = w[n - 1].pop(0)
                stack.append([a, b, c, 0, 0, 0])
                modQFA(a, b, c)
                w_tmp.append(c)

            w[n - 1].extend(w_tmp)

        for i in reversed(range(0, n - 1)):
            w_tmp = []
            w_tmp2 = []
            while (len(w[i]) > 2):
                a = w[i].pop(0)
                b = w[i].pop(0)
                c = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, c, carry])
                if TD == 2:
                    QFA(a, b, c, carry, and_ancilla[and_idx])
                else:
                    QFA(a, b, c, carry)
                and_idx += 1
                w_tmp.append(c)
                w_tmp2.append(carry)

            if flag == 1 and len(w[i]) == 2:
                a = w[i].pop(0)
                b = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, carry])
                if TD == 2:
                    QHA(a, b, carry, and_ancilla[and_idx])
                else:
                    QHA(a, b, carry)
                and_idx += 1
                w_tmp.append(b)
                w_tmp2.append(carry)

            w[i].extend(w_tmp)
            w[i + 1].extend(w_tmp2)

        and_idx = 0
        for i in range(len(w)):
            if (len(w[i]) > 3):
                t += 1
                flag = 1

    and_idx = 0
    t = 1
    while (t):
        n = len(w)
        t = 0
        w_tmp = []
        if (len(w[n - 1]) > 2):
            while (len(w[n - 1]) > 2):
                a = w[n - 1].pop(0)
                b = w[n - 1].pop(0)
                c = w[n - 1].pop(0)
                stack.append([a, b, c, 0, 0, 0])
                modQFA(a, b, c)
                w_tmp.append(c)

            w[n - 1].extend(w_tmp)

        if (len(w[n - 1]) == 2):
            return w, stack

        for i in reversed(range(0, n - 1)):
            w_tmp = []
            w_tmp2 = []
            while (len(w[i]) > 2):
                a = w[i].pop(0)
                b = w[i].pop(0)
                c = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, c, carry])
                if TD == 2:
                    QFA(a, b, c, carry, and_ancilla[and_idx])
                else:
                    QFA(a, b, c, carry)
                and_idx += 1
                w_tmp.append(c)
                w_tmp2.append(carry)

            if len(w[i]) == 2:
                a = w[i].pop(0)
                b = w[i].pop(0)
                carry = ancilla[num]
                num+=1
                an.append(carry)
                stack.append([a, b, carry])
                if TD == 2:
                    QHA(a, b, carry, and_ancilla[and_idx])
                else:
                    QHA(a, b, carry)
                and_idx += 1
                w_tmp.append(b)
                w_tmp2.append(carry)

            w[i].extend(w_tmp)
            w[i + 1].extend(w_tmp2)

        and_idx = 0
        for i in range(len(w)):
            if (len(w[i]) > 2):
                t += 1
    return w, stack

def apply_final_cnot(eng, w, output=0):
    idx = -1
    for i in range(len(w)):
        if len(w[i]) == 1:
            if output != 0:
                ancilla = output[idx]
                idx -= 1
            else:
                ancilla = eng.allocate_qubit()
            CNOT | (w[i].pop(), ancilla)
            w[i].append(ancilla)

    return w

def transfer_to_rca_components(w):
    rca_A = []
    rca_B = []
    R = []
    while (len(w)):
        if len(w[0]) == 1:
            if (len(w) == 1):
                rca_A.append(w.pop(0).pop())
                rca_B.append(eng.allocate_qubit())
            else:
                R.append(w.pop(0).pop())

        elif len(w[0]) == 2:
            tmp = w.pop(0)
            rca_A.append(tmp.pop(0))
            rca_B.append(tmp.pop())

    return rca_A, rca_B, R

def reverse_adders(stack):
    while (len(stack)):
        tmp = stack.pop()
        if (len(tmp) == 2):
            CNOT | (tmp[0], tmp[1])
        if (len(tmp) == 3):
            QHA_R(tmp[0], tmp[1], tmp[2])
        if (len(tmp) == 4):
            QFA_R(tmp[0], tmp[1], tmp[2], tmp[3])
        if (len(tmp) == 5):
            CNOT | (tmp[0], tmp[1])
        if (len(tmp) == 6):
            modQFA_R(tmp[0], tmp[1], tmp[2])

def csa0(eng, input, input2, ancillas, outputs):
    w1 = initialize_circuit_and_qubits(input)
    w1, stack0, input2, num = apply_adders_and_handle_carry0(eng, w1, input2, ancillas[0][0], ancillas[0][1])
    w2 = initialize_circuit_and_qubits(input2)

    w1, stack1 = apply_adders_and_handle_carry1(eng, w1, stack0, ancillas[0][0], num, ancillas[0][1])
    w2, stack2 = apply_adders_and_handle_carry(eng, w2, ancillas[1][0], ancillas[1][1])

    w1 = apply_final_cnot(eng, w1, outputs[0])
    w2 = apply_final_cnot(eng, w2, outputs[1])

    rca_A1, rca_B1, R1 = transfer_to_rca_components(w1)
    rca_A2, rca_B2, R2 = transfer_to_rca_components(w2)

    result1 = outDraper(eng, rca_A1, rca_B1, ancillas[0][1], outputs[0])
    result2 = outDraper(eng, rca_A2, rca_B2, ancillas[1][1], outputs[1])

    result1 = R1+result1
    result2 = R2+result2

    reverse_adders(stack2)
    reverse_adders(stack1)

    return result1, result2

def csa0_0(eng, input, input2, ancillas, outputs):
    w1 = initialize_circuit_and_qubits(input)
    w1, stack0, input2, num = apply_adders_and_handle_carry00(eng, w1, input2, ancillas[0][0]+ancillas[1][0][:31], ancillas[0][1])
    w2 = initialize_circuit_and_qubits(input2)

    w1, stack1 = apply_adders_and_handle_carry01(eng, w1, stack0, ancillas[0][0], num, ancillas[0][1])
    w2, stack2 = apply_adders_and_handle_carry(eng, w2, ancillas[1][0][31:], ancillas[1][1])

    w1 = apply_final_cnot(eng, w1, outputs[0])
    w2 = apply_final_cnot(eng, w2, outputs[1])

    rca_A1, rca_B1, R1 = transfer_to_rca_components(w1)
    rca_A2, rca_B2, R2 = transfer_to_rca_components(w2)

    result1 = outDraper(eng, rca_A1, rca_B1, ancillas[0][1], outputs[0])
    result2 = outDraper(eng, rca_A2, rca_B2, ancillas[1][1], outputs[1])

    result1 = R1+result1
    result2 = R2+result2

    reverse_adders(stack2)
    reverse_adders(stack1)

    return result1, result2

def csa1(eng, input, input2, input3, ancillas, outputs):
    w1 = initialize_circuit_and_qubits(input)
    if len(input3) != 2:
        w3 = initialize_circuit_and_qubits(input3)
    w1, stack0, input2, num = apply_adders_and_handle_carry0(eng, w1, input2, ancillas[0][0], ancillas[0][1])
    w2 = initialize_circuit_and_qubits(input2)

    w1, stack1 = apply_adders_and_handle_carry1(eng, w1, stack0, ancillas[0][0], num, ancillas[0][1])
    w2, stack2 = apply_adders_and_handle_carry(eng, w2, ancillas[1][0], ancillas[1][1])
    if len(input3) != 2:
        w3, stack3 = apply_adders_and_handle_carry(eng, w3, ancillas[2][0], ancillas[2][1])

    w1 = apply_final_cnot(eng, w1, outputs[0])
    w2 = apply_final_cnot(eng, w2, outputs[1])
    if len(input3) != 2:
        w3 = apply_final_cnot(eng, w3, outputs[2])

    rca_A1, rca_B1, R1 = transfer_to_rca_components(w1)
    rca_A2, rca_B2, R2 = transfer_to_rca_components(w2)
    if len(input3) != 2:
        rca_A3, rca_B3, R3 = transfer_to_rca_components(w3)

    result1 = outDraper(eng, rca_A1, rca_B1, ancillas[0][1], outputs[0])
    result2 = outDraper(eng, rca_A2, rca_B2, ancillas[1][1], outputs[1])
    if len(input3) != 2:
        result3 = outDraper(eng, rca_A3, rca_B3, ancillas[2][1], outputs[2])
    else:
        result3 = outDraper(eng, input3[0], input3[1], ancillas[2][0]+ancillas[2][1], outputs[2])

    result1 = R1+result1
    result2 = R2+result2
    if len(input3) != 2:
        result3 = R3+result3

    if len(input3) != 2:
        reverse_adders(stack3)
    reverse_adders(stack2)
    reverse_adders(stack1)

    return result1, result2, result3

def csa1_0(eng, input, input2, input3, ancillas, outputs):
    w1 = initialize_circuit_and_qubits(input)
    if len(input3) != 2:
        w3 = initialize_circuit_and_qubits(input3)
    w1, stack0, input2, num = apply_adders_and_handle_carry00(eng, w1, input2, ancillas[0][0]+ancillas[1][0][:2], ancillas[0][1])
    w2 = initialize_circuit_and_qubits(input2)

    w1, stack1 = apply_adders_and_handle_carry01(eng, w1, stack0, ancillas[0][0], num, ancillas[0][1])
    w2, stack2 = apply_adders_and_handle_carry(eng, w2, ancillas[1][0][2:], ancillas[1][1])
    if len(input3) != 2:
        w3, stack3 = apply_adders_and_handle_carry(eng, w3, ancillas[2][0], ancillas[2][1])

    w1 = apply_final_cnot(eng, w1, outputs[0])
    w2 = apply_final_cnot(eng, w2, outputs[1])
    if len(input3) != 2:
        w3 = apply_final_cnot(eng, w3, outputs[2])

    rca_A1, rca_B1, R1 = transfer_to_rca_components(w1)
    rca_A2, rca_B2, R2 = transfer_to_rca_components(w2)
    if len(input3) != 2:
        rca_A3, rca_B3, R3 = transfer_to_rca_components(w3)

    result1 = outDraper(eng, rca_A1, rca_B1, ancillas[0][1], outputs[0])
    result2 = outDraper(eng, rca_A2, rca_B2, ancillas[1][1], outputs[1])
    if len(input3) != 2:
        result3 = outDraper(eng, rca_A3, rca_B3, ancillas[2][1], outputs[2])
    else:
        result3 = outDraper(eng, input3[0], input3[1], ancillas[2][0]+ancillas[2][1], outputs[2])
    result1 = R1+result1
    result2 = R2+result2
    if len(input3) != 2:
        result3 = R3+result3

    if len(input3) != 2:
        reverse_adders(stack3)
    reverse_adders(stack2)
    reverse_adders(stack1)

    return result1, result2, result3

def csa1_d(eng, input, input2, input3, ancillas, outputs):
    w1 = initialize_circuit_and_qubits(input)
    if len(input3) != 2:
        w3 = initialize_circuit_and_qubits(input3)
    w1, stack0, input2, num = apply_adders_and_handle_carry0(eng, w1, input2, ancillas[0][0], ancillas[0][1])
    w2 = initialize_circuit_and_qubits(input2)

    w1, stack1 = apply_adders_and_handle_carry1(eng, w1, stack0, ancillas[0][0], num, ancillas[0][1])
    w2, stack2 = apply_adders_and_handle_carry(eng, w2, ancillas[1][0], ancillas[1][1])
    if len(input3) != 2:
        w3, stack3 = apply_adders_and_handle_carry(eng, w3, ancillas[2][0], ancillas[2][1])

    w1 = apply_final_cnot(eng, w1, outputs[0])
    w2 = apply_final_cnot(eng, w2, outputs[1])
    if len(input3) != 2:
        w3 = apply_final_cnot(eng, w3, outputs[2])

    rca_A1, rca_B1, R1 = transfer_to_rca_components(w1)
    rca_A2, rca_B2, R2 = transfer_to_rca_components(w2)
    if len(input3) != 2:
        rca_A3, rca_B3, R3 = transfer_to_rca_components(w3)

    result1 = outDraper(eng, rca_A1, rca_B1, ancillas[0][1], outputs[0])
    result2 = outDraper(eng, rca_A2, rca_B2, ancillas[1][1], outputs[1])
    if len(input3) != 2:
        result3 = outDraper_dag(eng, rca_A3, rca_B3, ancillas[2][1], outputs[2])
    else:
        result3 = outDraper_dag(eng, input3[0], input3[1], ancillas[2][0]+ancillas[2][1], outputs[2])

    result1 = R1+result1
    result2 = R2+result2
    if len(input3) != 2:
        result3 = R3+result3

    if len(input3) != 2:
        reverse_adders(stack3)
    reverse_adders(stack2)
    reverse_adders(stack1)

    return result1, result2, result3

def csa2(eng, input, ancillas, output):
    w = initialize_circuit_and_qubits(input)
    w, stack = apply_adders_and_handle_carry(eng, w, ancillas[0][0], ancillas[0][1])
    w = apply_final_cnot(eng, w, output)
    rca_A, rca_B, R = transfer_to_rca_components(w)
    result = outDraper(eng, rca_A, rca_B, ancillas[0][1], output)
    reverse_adders(stack)

    return R+result

def csa2_d(eng, input, ancillas, output): # dag
    if len(input) != 2:
        w = initialize_circuit_and_qubits(input)
        w, stack = apply_adders_and_handle_carry(eng, w, ancillas[0][0], ancillas[0][1])
        w = apply_final_cnot(eng, w, output)
        rca_A, rca_B, R = transfer_to_rca_components(w)
        result = outDraper_dag(eng, rca_A, rca_B, ancillas[0][1], output)
        reverse_adders(stack)
        result = R+result
    else:
        result = outDraper_dag(eng, input[0], input[1], ancillas[0][0]+ancillas[0][1], output)

    return result

def make_list(num,length):

    arr1 = list(format(num, 'b'))
    if length > len(arr1):
        arr2 = ['0' for i in range(length-len(arr1))]
        arr2.extend(arr1)
        arr1 = arr2
    return arr1

def SHA256_init_H(H, idx):

    H0_hex = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

    chain = make_list(H0_hex[idx], 32)
    for j in range(32):
        if chain[j] == '1':
            # X | H[i][j]
            X | H[idx][31-j]

    return H[idx]

def SHA256_init_K(K, idx):

    K_hex = [0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
             0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
             0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
             0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
             0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
             0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
             0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
             0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
             0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
             0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
             0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2]

    chain = make_list(K_hex[idx], 32)
    for j in range(32):
        if chain[j] == '1':
            # X | K[i][j]
            X | K[31-j]


def SHA256_preprocessing(msg, q_msg):
    global msg_bin_arr

    # Padding
    lenght = len(msg)
    msg_bin = []

    for i in range(lenght):
        msg_bin += make_list(ord(msg[i]), 8)
    lenght = len(msg_bin)
    msg_bin += ['1']

    for i in range(lenght+1, 448):
        msg_bin += ['0']

    len_arr = make_list(lenght, 64)
    msg_bin += len_arr

    msg_bin_arr = msg_bin

    if resource_check != 1:
        for i in range(16):
            init_Ws(q_msg, i)

def init_W(W, idx):
    # Parsing
    for i in range(32):
        if msg_bin_arr[idx * 32 + i] == '1':
            X | W[31-i]

def init_Ws(W, idx):
    # Parsing
    for i in range(32):
        if msg_bin_arr[idx * 32 + i] == '1':
            X | W[idx][31-i]

def w_sigma(rr0, rr1, ss, input, output):

    ss_len = 32 - ss
    for i in range(32):
        CNOT | (input[(rr0+i)%32], output[i])
    for i in range(32):
        CNOT | (input[(rr1+i)%32], output[i])
    for i in range(32):
        if i < ss_len: # ss ~ 31
            CNOT | (input[i+ss], output[i])

def r_sigma(rr0, rr1, rr2, input, output):

    for i in range(32):
        CNOT | (input[(rr0+i)%32], output[i])
    for i in range(32):
        CNOT | (input[(rr1+i)%32], output[i])
    for i in range(32):
        CNOT | (input[(rr2+i)%32], output[i])

def maj(x, y, z, tmp = 0, ancilla=0):
    for i in range(32):
        CNOT | (z[i], y[i])
    for i in range(32):
        CNOT | (z[i], x[i])
    if TD == 2:
        for i in range(32):
            toffoli_gate(x[i], y[i], tmp[i], ancilla[i])
            CNOT | (tmp[i], z[i])
    else:
        for i in range(32):
            toffoli_gate(x[i], y[i], z[i])

def maj_dag(x, y, z, tmp = 0):
    if TD == 2:
        for i in range(32):
            CNOT | (tmp[i], z[i])
            toffoli_gate(x[i], y[i], tmp[i], 0, False)
    else:
        for i in range(32):
            toffoli_gate(x[i], y[i], z[i])
    for i in range(32):
        CNOT | (z[i], x[i])
    for i in range(32):
        CNOT | (z[i], y[i])

def ch(x, y, z, tmp=0, ancilla=0):
    for i in range(32):
        CNOT | (z[i], y[i])
    if TD == 2:
        for i in range(32):
            toffoli_gate(x[i], y[i], tmp[i], ancilla[i])
            CNOT | (tmp[i], z[i])
    else:
        for i in range(32):
            toffoli_gate(x[i], y[i], z[i])

def ch_dag(x, y, z, tmp=0):
    if TD == 2:
        for i in range(32):
            CNOT | (tmp[i], z[i])
            toffoli_gate(x[i], y[i], tmp[i], 0, False)
    else:
        for i in range(32):
            toffoli_gate(x[i], y[i], z[i])
    for i in range(32):
        CNOT | (z[i], y[i])

def copy(x, y):
    for i in range(32):
        CNOT | (x[i], y[i])

def SHA256_round(K, W, H, ancilla, add_ancilla, garbages):
    ancilla_0 = 0
    ancilla_1 = 0

    gg1 = 0
    gg2 = 0

    h = []
    for i in range(8):
        h.append(H[i])

    gar_idx_d = 27
    dag_idx = 43
    init_idx = 8
    idx = 16  # W index

    ## round proceess ##
    for i in range(64):
        if i == 0:
            w_sigma(7, 18, 3, W[idx - 15], garbages[-7])  # sigma0
        elif i >= 1 and i <= 39:
            if i <= 7 or i >= 14:
                w_sigma(7, 18, 3, W[idx - 15], ancilla_0)  # sigma0
            w_sigma(17, 19, 10, W[idx - 2], ancilla_1)  # sigma1
            if i >= 32:
                w_sigma(7, 18, 3, W[(idx+1) - 15], garbages[-7])  # sigma0
                w_sigma(17, 19, 10, W[(idx+1) - 2], garbages[-6])  # sigma1
        elif i >= 44 and i <= 49:
            w_sigma(7, 18, 3, W[dag_idx - 15], ancilla_0)  # sigma0
            w_sigma(17, 19, 10, W[dag_idx - 2], ancilla_1)  # sigma1
            w_sigma(7, 18, 3, W[(dag_idx-1) - 15], garbages[-7])  # sigma0
            w_sigma(17, 19, 10, W[(dag_idx-1) - 2], garbages[-6])  # sigma1
        elif i >= 50 and i <= 57:
            if i == 50: ## -8:15
                init_W(garbages[-8], dag_idx-16) # W[15]
                w_sigma(7, 18, 3, W[dag_idx-15], ancilla_0)  # sigma0
                w_sigma(17, 19, 10, W[dag_idx-2], ancilla_1)  # sigma1
                w_sigma(7, 18, 3, garbages[-8], garbages[-7])  # sigma0
                w_sigma(17, 19, 10, W[(dag_idx-1)-2], garbages[-6])  # sigma1
            elif i >= 51 and i <= 53:
                if i == 53: ## -8:8
                    init_W(garbages[-8], dag_idx - 16 -1)  # W[8]
                    w_sigma(7, 18, 3, garbages[-8], ancilla_0)
                w_sigma(17, 19, 10, W[dag_idx - 2], ancilla_1)  # sigma1
                w_sigma(17, 19, 10, W[(dag_idx - 1)-2], garbages[-6])  # sigma1
            elif i == 54:
                init_W(garbages[-8], dag_idx - 1 - 7)  # W[15]
                w_sigma(17, 19, 10, W[dag_idx - 2], ancilla_1)  # sigma1
                w_sigma(7, 18, 3, W[7], garbages[-7])  # sigma0
                w_sigma(17, 19, 10, W[(dag_idx - 1)-2], garbages[-6])  # sigma1
            elif i == 55 or i == 56:
                if i == 55:
                    init_W(garbages[-8], dag_idx - 15 + 2)  # W[8] dag
                w_sigma(7, 18, 3, W[dag_idx-15], ancilla_0) # sigma0
                w_sigma(17, 19, 10, W[dag_idx - 2], ancilla_1)  # sigma1
                w_sigma(7, 18, 3, W[(dag_idx-15)-1], garbages[-7])  # sigma0
                w_sigma(17, 19, 10, W[(dag_idx - 1) - 2], garbages[-6])  # sigma1
            elif i == 57: ## -8:15
                init_W(garbages[-8], dag_idx -2)  # W[15]
                w_sigma(17, 19, 10, garbages[-8], ancilla_1)  # sigma1
                w_sigma(7, 18, 3, W[2], ancilla_0)  # sigma0
                w_sigma(7, 18, 3, W[1], garbages[-7])  # sigma0
        elif i == 58:
            final = garbages[-10:]
            idx = 0

        SHA256_init_K(K,i)
        r_sigma(2, 13, 22, h[0], ancilla[0]) # 32
        r_sigma(6, 11, 25, h[4], ancilla[1]) # 32

        if i <= 57:
            maj(h[0], h[1], h[2], garbages[-9], add_ancilla[1][0])
            ch(h[4],h[5],h[6], garbages[-10], add_ancilla[1][1])
        else:
            maj(h[0], h[1], h[2], ancilla_0, add_ancilla[1][0])
            ch(h[4], h[5], h[6], ancilla_1, add_ancilla[1][1])

        if i >= 9 and i <= 14:
            input = [ancilla[1], h[6], h[7], K, ancilla[0], h[2]]
            input2 = [K, h[3]]
        else:
            input = [ancilla[1], h[6], h[7], W[i], K, ancilla[0], h[2]]
            input2 = [W[i], K, h[3]]

        tmp1 = garbages[-1] + garbages[-2]
        tmp2 = garbages[-3] + garbages[-4] + garbages[-5]
        tmp = [[tmp1, tmp2]]

        if i == 0:
            input3 = [W[idx - 16], garbages[-7]]
            output = [0,0,garbages[0]]
            h[7], h[3], ww = csa1(eng, input, input2, input3, add_ancilla, output)
            W.append(ww)

        elif i >= 1 and i <= 31:
            if i >= 1 and i <= 5:
                input3 = [W[idx - 16], ancilla_0, ancilla_1]
            elif i == 6 or i == 7:
                input3 = [W[idx - 16], ancilla_0, W[idx - 7], ancilla_1]
            elif i == 8:
                input3 = [W[idx - 16], W[idx - 7], ancilla_1]
            elif i >= 9 and i <= 13:
                input3 = [W[idx - 7], ancilla_1]
            elif i == 14:
                input3 = [ancilla_0, W[idx - 7], ancilla_1]
            else:
                input3 = [W[idx - 16], ancilla_0, W[idx - 7], ancilla_1]

            if i <= 27:
                output = [0, 0, garbages[i]]
            else:
                output = [0, 0, 0]

            if i >= 9 and i <= 14:
                h[7], h[3], ww = csa1_0(eng, input, input2, input3, add_ancilla, output)
            else:
                h[7], h[3], ww = csa1(eng, input, input2, input3, add_ancilla, output)
            W.append(ww)

            if i >= 24:
                init_Ws(W, i - 16)

        elif i >= 32 and i <= 39:
            input3 = [W[idx - 16], ancilla_0, W[idx - 7], ancilla_1]
            if i <= 35:
                output = [0, 0, 0]
                input4 = [W[(idx+1)-16], garbages[-7], W[(idx+1)-7], garbages[-6]]
                h[7], h[3], ww1 = csa1(eng, input, input2, input3, add_ancilla, output)
                ww2 = csa2(eng, input4, tmp, 0)
                W.append(ww1)
                W.append(ww2)

            else:
                output = [0, 0, W[init_idx]]
                input4 = [W[(idx+1)-16], garbages[-7], W[(idx+1)-7], garbages[-6]]
                h[7], h[3], ww1 = csa1(eng, input, input2, input3, add_ancilla, output)
                ww2 = csa2(eng, input4, tmp, W[init_idx+1])
                W.append(ww1)
                W.append(ww2)
                init_idx += 2

        elif i >= 40 and i <= 43:
            output = [0, 0]
            h[7], h[3] = csa0(eng, input, input2, add_ancilla[0:2], output)

        elif i >= 44 and i <= 49: # gar ancilla dagger
            if i == 44:
                input3 = [W[dag_idx - 16], ancilla_0, W[dag_idx - 7], ancilla_1]
                output = [0, 0, garbages[gar_idx_d]]
                input4 = [W[(dag_idx-1) - 16], garbages[-7], W[(dag_idx-1) - 7], garbages[-6]]
                gg2 = csa2_d(eng, input4, tmp, garbages[gar_idx_d-1])
                h[7], h[3], gg1 = csa1_d(eng, input, input2, input3, add_ancilla, output)
            else:
                input3 = [W[dag_idx - 16], ancilla_0, W[dag_idx - 7], ancilla_1]
                output = [gg1, gg2, garbages[gar_idx_d]]
                input4 = [W[(dag_idx-1) - 16], garbages[-7], W[(dag_idx-1) - 7], garbages[-6]]
                h[7], h[3], gg1 = csa1_d(eng, input, input2, input3, add_ancilla, output)
                gg2 = csa2_d(eng, input4, tmp, garbages[gar_idx_d-1])
            gar_idx_d -= 2

        elif i >= 50 and i <= 57: # gar ancilla dagger
            if i == 50:
                input3 = [garbages[-8], ancilla_0, W[dag_idx-7], ancilla_1]
                output = [gg1, gg2, garbages[gar_idx_d]]
                input4 = [garbages[-7], W[(dag_idx-1)-7], garbages[-6]]
                h[7], h[3], gg1 = csa1_d(eng, input, input2, input3, add_ancilla, output)
                gg2 = csa2_d(eng, input4, tmp, garbages[gar_idx_d - 1])
                gar_idx_d -= 2
            elif i >= 51 and i <= 53:
                if i == 53:
                    input4 = [garbages[-8], W[(dag_idx-1)-7], garbages[-6]]
                else:
                    input4 = [W[(dag_idx-1)-7], garbages[-6]]
                input3 = [W[dag_idx - 7], ancilla_1]
                output = [gg1, gg2, garbages[gar_idx_d]]
                h[7], h[3], gg1 = csa1_d(eng, input, input2, input3, add_ancilla, output)
                gg2 = csa2_d(eng, input4, tmp, garbages[gar_idx_d - 1])
                gar_idx_d -= 2
            elif i == 54: # -8:15
                input3 = [W[7], ancilla_0, W[dag_idx-7], ancilla_1]
                output = [gg1, gg2, garbages[gar_idx_d]]
                input4 = [W[6], garbages[-7], garbages[-8], garbages[-6]]
                h[7], h[3], gg1 = csa1_d(eng, input, input2, input3, add_ancilla, output)
                gg2 = csa2_d(eng, input4, tmp, garbages[gar_idx_d - 1])
                gar_idx_d -= 2
            elif i == 55 or i == 56:
                input3 = [W[dag_idx-16], ancilla_0, ancilla_1]
                output = [gg1, gg2, garbages[gar_idx_d]]
                input4 = [W[(dag_idx-16)-1], garbages[-7], garbages[-6]]
                h[7], h[3], gg1 = csa1_d(eng, input, input2, input3, add_ancilla, output)
                gg2 = csa2_d(eng, input4, tmp, garbages[gar_idx_d - 1])
                gar_idx_d -= 2
            elif i == 57: ## -8:15
                input3 = [W[1], ancilla_0, ancilla_1]
                output = [gg1, gg2, garbages[gar_idx_d]]
                input4 = [W[0], garbages[-7]]
                h[7], h[3], gg1 = csa1_d(eng, input, input2, input3, add_ancilla, output)
                gg2 = csa2_d(eng, input4, tmp, garbages[gar_idx_d - 1])
                gar_idx_d -= 2
        elif i >= 58:
            if i == 58:
                output = [gg1, gg2]
            else:
                output = [final[idx], final[idx+1]]
                idx += 2

            h[7], h[3] = csa0(eng, input, input2, add_ancilla[0:2], output)

        if i == 0:
            ancilla_0 = SHA256_init_H(H, 3)  # dag
            ancilla_1 = SHA256_init_H(H, 7)  # dag
            w_sigma(7, 18, 3, W[idx - 15], garbages[-7])  # dag
            idx += 1

        elif i >= 1 and i <= 39:
            w_sigma(17, 19, 10, W[idx - 2], ancilla_1)  # sigma1
            if i <= 7 or i >= 14:
                w_sigma(7, 18, 3, W[idx - 15], ancilla_0)  # sigma0
            if i >= 32:
                w_sigma(17, 19, 10, W[(idx+1) - 2], garbages[-6])  # sigma1 dag
                w_sigma(7, 18, 3, W[(idx+1) - 15], garbages[-7])  # sigma0 dag
                idx += 1
            idx += 1

        elif i >= 44 and i <= 49:
            w_sigma(17, 19, 10, W[(dag_idx-1) - 2], garbages[-6])  # sigma1 dag
            w_sigma(7, 18, 3, W[(dag_idx-1) - 15], garbages[-7])  # sigma0 dag
            w_sigma(17, 19, 10, W[dag_idx - 2], ancilla_1)  # sigma1 dag
            w_sigma(7, 18, 3, W[dag_idx - 15], ancilla_0)  # sigma0 dag
            dag_idx -= 2

        elif i >= 50 and i <= 57:
            if i == 50:
                w_sigma(17, 19, 10, W[(dag_idx-1)-2], garbages[-6])  # sigma1
                w_sigma(7, 18, 3, garbages[-8], garbages[-7])  # sigma0
                w_sigma(17, 19, 10, W[dag_idx-2], ancilla_1)  # sigma1
                w_sigma(7, 18, 3, W[dag_idx-15], ancilla_0)  # sigma0
                init_W(garbages[-8], dag_idx-16) # W[15] dag
                dag_idx -= 2
            elif i >= 51 and i <= 53:
                w_sigma(17, 19, 10, W[(dag_idx - 1) - 2], garbages[-6])  # sigma1
                w_sigma(17, 19, 10, W[dag_idx - 2], ancilla_1)  # sigma1
                if i == 53:
                    init_W(garbages[-8], dag_idx-16-1) # W[8] dag
                dag_idx -= 2
            elif i == 54:
                init_W(garbages[-8], dag_idx - 1 - 7)  # W[15] dag
                w_sigma(17, 19, 10, W[(dag_idx - 1) - 2], garbages[-6])  # sigma1
                w_sigma(7, 18, 3, W[7], garbages[-7])  # sigma0
                w_sigma(17, 19, 10, W[dag_idx - 2], ancilla_1)  # sigma1
                init_W(garbages[-8], dag_idx - 1 - 7 - 7)  # W[8] init
                w_sigma(7, 18, 3, garbages[-8], ancilla_0)  # sigma0 dag
                dag_idx -= 2
            elif i == 55 or i == 56:
                w_sigma(17, 19, 10, W[(dag_idx - 1) - 2], garbages[-6])  # sigma1
                w_sigma(7, 18, 3, W[(dag_idx - 15) - 1], garbages[-7])  # sigma0
                w_sigma(17, 19, 10, W[dag_idx - 2], ancilla_1)  # sigma1
                w_sigma(7, 18, 3, W[dag_idx - 15], ancilla_0)  # sigma0
                dag_idx -= 2
            elif i == 57:
                w_sigma(7, 18, 3, W[1], garbages[-7])  # sigma0
                w_sigma(7, 18, 3, W[2], ancilla_0)  # sigma0
                w_sigma(17, 19, 10, garbages[-8], ancilla_1)  # sigma1
                init_W(garbages[-8], dag_idx - 2)  # W[15] dag
                dag_idx -= 2

        if i <= 57:
            ch_dag(h[4],h[5],h[6], garbages[-10])
            maj_dag(h[0], h[1], h[2], garbages[-9])
        else:
            ch_dag(h[4], h[5], h[6], ancilla_1)
            maj_dag(h[0], h[1], h[2], ancilla_0)
        r_sigma(2, 13, 22, h[0], ancilla[0]) # dag
        r_sigma(6, 11, 25, h[4], ancilla[1]) # dag
        SHA256_init_K(K,i) # dag

        swap = h[0]
        h[0] = h[7]
        h[7] = h[6]
        h[6] = h[5]
        h[5] = h[4]
        h[4] = h[3]
        h[3] = h[2]
        h[2] = h[1]
        h[1] = swap

        if resource_check != 1:
            print('\n')
            print("round ", i)
            for k in range(8):
                print_vector(h[k], 32)
        else:
            print("round ", i)

def print_vector(element, length):

    All(Measure) | element
    for k in range(length):
        print(int(element[length-1-k]), end='')
    print()

def sha2(eng):
    H = []
    K = eng.allocate_qureg(32)
    W = []

    for i in range(8):
        H.append(eng.allocate_qureg(32))

    for i in range(16):
        W.append(eng.allocate_qureg(32))

    ancilla = []
    for i in range(2):
        ancilla.append(eng.allocate_qureg(32))  # r_sigma

    add_ancilla = []
    tmp = []
    tmp.append(eng.allocate_qureg(153)) # 7 CSA
    tmp.append(eng.allocate_qureg(75)) # outDraper_30bit
    add_ancilla.append(tmp)
    tmp = []
    tmp.append(eng.allocate_qureg(92)) # 5 CSA
    tmp.append(eng.allocate_qureg(75)) # outDraper_30bit
    add_ancilla.append(tmp)
    tmp = []
    tmp.append(eng.allocate_qureg(62))  # 5 CSA
    tmp.append(eng.allocate_qureg(78))  # outDraper_31bit
    add_ancilla.append(tmp)

    garbages = []
    for i in range(38):
        garbages.append(eng.allocate_qureg(32))

    msg = "abc"

    for i in range(8):
        SHA256_init_H(H, i)
    SHA256_preprocessing(msg, W)
    SHA256_round(K, W, H, ancilla, add_ancilla, garbages)
    ### Constant addition operation performed at the final stage is omitted ###


global TD
global resource_check

TD = 0
resource_check = 0
print('Generate Ciphertext...')
Simulate = ClassicalSimulator()
eng = MainEngine(Simulate)
sha2(eng)
print('\n')
eng.flush()

TD = 1
resource_check = 1
print('Estimate cost...')
Resource = ResourceCounter()
eng = MainEngine(Resource)
sha2(eng)
print(Resource)
print('\n')
eng.flush()


## AND gate ##
TD = 2
resource_check = 1
print('[AND] Estimate cost...')
Resource = ResourceCounter()
eng = MainEngine(Resource)
sha2(eng)
print(Resource)
print('\n')
eng.flush()