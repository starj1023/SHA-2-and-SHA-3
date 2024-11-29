from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control, Dagger

def Round_constant_XOR(eng, k, rc, bit):
    print(rc)
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def print_state(eng, state):
    for i in range(25):
        for j in range(64):
            Measure | state[i][63-j]
            print((int(state[i][63 - j])), end='')
            if ((j + 1) % 4 == 0):
                print(" ", end='')
        print(" ")

def print_bc(eng, bc):
    for i in range(5):
        for j in range(64):
            Measure | bc[4-i][63-j]
            print((int(bc[4-i][63 - j])), end='')
            if ((j + 1) % 4 == 0):
                print(" ", end='')
        print(" ")

def print_input(eng, b, k):
    All(Measure) | b
    All(Measure) | k
    print('Plaintext : 0x', end='')
    print_hex(eng, b)
    print('\nKey : 0x', end='')
    print_hex(eng, k)
    print('\n')

def print_hex(eng, qubits):

    for i in reversed(range(32)):
        temp = 0
        temp = temp+int(qubits[4*i+3])*8
        temp = temp+int(qubits[4*i+2])*4
        temp = temp+int(qubits[4*i+1])*2
        temp = temp+int(qubits[4*i])

        temp = hex(temp)
        y = temp.replace("0x", "")
        print(y, end='')

def AND_gate(eng, a, b, c, ancilla):
    if (resource_check):
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
    else:
        Toffoli | (a, b, c)

def AND_gate_dag(eng, a, b, c):
    H | b
    CNOT | (a, b)
    X | c

    with Dagger(eng):
        Measure | c
    H | b
    H | c


def ROL(eng, state, keccakf_rotc):
    # keccakf_rotc = [1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44]
    out = []
    for i in range(keccakf_rotc):
        out.append(state[64-keccakf_rotc+i])

    for i in range(64 - keccakf_rotc):
        out.append(state[i])

    return out

def SHA(eng):

    keccakf_rotc = [1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44]
    keccakf_piln = [10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4, 15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1]
    keccakf_rndc = [0x0000000000000001, 0x0000000000008082, 0x800000000000808a,
        0x8000000080008000, 0x000000000000808b, 0x0000000080000001,
        0x8000000080008081, 0x8000000000008009, 0x000000000000008a,
        0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
        0x000000008000808b, 0x800000000000008b, 0x8000000000008089,
        0x8000000000008003, 0x8000000000008002, 0x8000000000000080,
        0x000000000000800a, 0x800000008000000a, 0x8000000080008081,
        0x8000000000008080, 0x0000000080000001, 0x8000000080008008]

    rows = 25
    cols = 64
    state = [[0 for j in range(cols)] for i in range(rows)]
    for i in range(rows):
        for j in range(cols):
            state[i][j] = eng.allocate_qubit()

    if (resource_check != 1):
        for i in range(25):
            for j in range(8):
                X | state[i][0 + 8 * j]
                X | state[i][1 + 8 * j]

                X | state[i][5 + 8 * j]
                X | state[i][7 + 8 * j]


    copy_state1 = [[0 for j in range(cols)] for i in range(rows)]
    for i in range(rows):
        for j in range(cols):
            copy_state1[i][j] = eng.allocate_qubit()

    copy_state1_second = [[0 for j in range(cols)] for i in range(rows)]
    for i in range(rows):
        for j in range(cols):
            copy_state1_second[i][j] = eng.allocate_qubit()

    for i in range(25):
        for j in range(64):
            X | copy_state1[i][j]

    for i in range(25):
        for j in range(64):
            X | copy_state1_second[i][j]

    kb_a = 5
    kb_b = 5
    kb_c = 64

    anc = [[[0 for k in range(kb_c)] for j in range(kb_b)] for i in range(kb_a)]
    for i in range(kb_a):
        for j in range(kb_b):
            for k in range(kb_c):
                anc[i][j][k] = eng.allocate_qubit()

    copy_state2 = [[[0 for _ in range(24)] for _ in range(25)] for _ in range(64)]

    bc = []
    bc = [[[0 for _ in range(24)] for _ in range(5)] for _ in range(64)]

    copy_state2_before_phi = [[[0 for _ in range(24)] for _ in range(25)] for _ in range(64)]

    bc[0] = generate_bc(eng)
    bc[1] = generate_bc(eng)
    bc[2] = generate_bc(eng)
    bc[3] = generate_bc(eng)
    bc[4] = generate_bc(eng)

    copy_state2[0] = generate_state(eng)
    copy_state2[1] = generate_state(eng)
    copy_state2[2] = generate_state(eng)
    copy_state2[3] = generate_state(eng)
    copy_state2[4] = generate_state(eng)

    copy_state2[5] = generate_state(eng)
    copy_state2[9] = generate_state(eng)
    copy_state2[13] = generate_state(eng)
    copy_state2[17] = generate_state(eng)
    copy_state2[21] = generate_state(eng)

    # 0 ~ 2 rounds
    for round in range(24):
        print('Round ', round)
        if (round != 0): #reverse
            state_copy(eng, copy_state2[round-1], target)

        if (round < 4 or (round > 7 and round < 12) or (round > 15 and round < 20)):
            reverse_target = copy_state1_second
        else:
            reverse_target = copy_state1

        index = int(round / 4) - 1
        jump = 4 * index

        # Clean ancilla qubits in Round 3, 2, 1, 0
        if (round >= 4 and (round-4*index) % 4 == 0):
            print('reverse3')
            reverse_theta(eng, bc[3+jump], copy_state2_before_phi[3+jump]) # bc[round-1] -> clean (C)
            bc[5+jump] = bc[3+jump]

        if (round >= 4 and (round-4*index) % 5 == 0):
            print('reverse2')
            state_copy(eng, copy_state2[2+jump], reverse_target)
            Round_constant_XOR(eng, copy_state2_before_phi[3+jump][0], keccakf_rndc[2+jump], 64)

            state_copy(eng, copy_state2[2 + jump],

            copy_state2_before_phi[3 + jump])  # copy_state2_before_phi[round - 1] -> clean (B)
            chi_reverse(eng, reverse_target, copy_state2[2+jump], copy_state2_before_phi[3+jump], anc)

            state_copy(eng, copy_state2[2+jump], reverse_target)

            copy_state2[6+jump] = copy_state2_before_phi[3+jump]

            reverse_theta(eng, bc[2+jump], copy_state2_before_phi[2+jump]) # bc[round-2] -> clean (B)
            bc[6+jump] = bc[2+jump]

        if (round >= 4 and (round-4*index) % 6 == 0):
            state_copy(eng, copy_state2[1+jump], reverse_target)

            Round_constant_XOR(eng, copy_state2_before_phi[2+jump][0], keccakf_rndc[1+jump], 64)
            state_copy(eng, copy_state2[1 + jump],
            copy_state2_before_phi[2 + jump])  # copy_state2_before_phi[round - 2] -> clean (B)
            chi_reverse(eng, reverse_target, copy_state2[1+jump], copy_state2_before_phi[2+jump], anc)

            state_copy(eng, copy_state2[1+jump], reverse_target)  # Clean copy_state (save 1600 qubits)

            copy_state2[7+jump] = copy_state2_before_phi[2+jump]

            reverse_theta(eng, bc[1+jump], copy_state2_before_phi[1+jump]) # bc[round-3] -> clean (B)
            bc[7+jump] = bc[1+jump]

        if (round >= 4 and (round-4*index) % 7 == 0):
            if(round==23):
                print('Skip reverse')
            else:
                print('Reverse')

                state_copy(eng, copy_state2[0+jump], reverse_target)

                Round_constant_XOR(eng, copy_state2_before_phi[1+jump][0], keccakf_rndc[0+jump], 64)
                state_copy(eng, copy_state2[0 + jump],
                copy_state2_before_phi[1 + jump])  # copy_state2_before_phi[round - 3] -> clean (A)
                chi_reverse(eng, reverse_target, copy_state2[0+jump], copy_state2_before_phi[1+jump], anc)

                state_copy(eng, copy_state2[0+jump], reverse_target)  # copy_state1 -> clean (A), note that 0b111111...

                copy_state2[8+jump] = copy_state2_before_phi[1+jump]
                #copy_state 1 is free in round 8

                reverse_theta(eng, bc[0+jump], copy_state2_before_phi[0+jump])  # bc[round-4] -> clean (A)
                bc[8+jump] = bc[0+jump]

        theta(eng, bc[round], state)

        copy_state2_before_phi[round] = store_before(eng, state)
        phi(eng, state, keccakf_piln, keccakf_rotc)

        if (round < 4 or (round > 7 and round < 12) or (round > 15 and round < 20)):
            target = copy_state1
        else:
            target = copy_state1_second

        state_copy(eng, state, target)

        chi(eng, target, state, copy_state2[round])

        state_copy(eng, state, copy_state2[round])

        temp = state
        state = copy_state2[round]
        copy_state2[round] = temp

        Round_constant_XOR(eng, state[0], keccakf_rndc[round], 64)

        if (resource_check != 1):
            print('Round state: ', round)
            print_state(eng, state)

def generate_bc(eng):
    rows = 5
    cols = 64
    bc = [[0 for j in range(cols)] for i in range(rows)]
    for i in range(rows):
        for j in range(cols):
            bc[i][j] = eng.allocate_qubit()

    return bc

def store_before(eng, state):
    rows = 25
    cols = 64
    before = [[0 for j in range(cols)] for i in range(rows)]
    for i in range(rows):
        for j in range(cols):
            before[i][j] = state[i][j]

    return before

def generate_state(eng):
    rows = 25
    cols = 64
    state = [[0 for j in range(cols)] for i in range(rows)]
    for i in range(rows):
        for j in range(cols):
            state[i][j] = eng.allocate_qubit()

    return state

def state_copy_x(eng, state, state_copy):
    for i in range(25):
        for j in range(64):
            CNOT | (state[i][j], state_copy[i][j])
            X | state_copy[i][j]

def state_copy(eng, state, state_copy):
    for i in range(25):
        for j in range(64):
            CNOT | (state[i][j], state_copy[i][j])

def theta(eng, bc, state):
    for i in range(5):
        for j in range(64):
            CNOT | (state[i][j], bc[i][j])
            CNOT | (state[i + 5][j], bc[i][j])
            CNOT | (state[i + 10][j], bc[i][j])
            CNOT | (state[i + 15][j], bc[i][j])
            CNOT | (state[i + 20][j], bc[i][j])

    j = 0
    for i in range(5):
        for k in range(64):
            CNOT | (bc[(i + 4) % 5][k], state[j * 5 + i][k])
            CNOT | (bc[(i + 1) % 5][(k - 1) % 64], state[(j + 1) * 5 + i][k])
            CNOT | (bc[(i + 1) % 5][(k - 1) % 64], state[(j + 2) * 5 + i][k])
            CNOT | (bc[(i + 1) % 5][(k - 1) % 64], state[(j + 3) * 5 + i][k])
            CNOT | (bc[(i + 1) % 5][(k - 1) % 64], state[(j + 4) * 5 + i][k])

    for i in range(5):
        for k in range(64):
            CNOT | (bc[(i + 1) % 5][(k - 1) % 64], state[j * 5 + i][k])
            CNOT | (bc[(i + 4) % 5][k], state[(j + 1) * 5 + i][k])
            CNOT | (bc[(i + 4) % 5][k], state[(j + 2) * 5 + i][k])
            CNOT | (bc[(i + 4) % 5][k], state[(j + 3) * 5 + i][k])
            CNOT | (bc[(i + 4) % 5][k], state[(j + 4) * 5 + i][k])

def reverse_theta(eng, bc, state):
    j = 0
    for i in range(5):
        for k in range(64):
            CNOT | (bc[(i + 4) % 5][k], state[j * 5 + i][k])
            CNOT | (bc[(i + 1) % 5][(k - 1) % 64], state[(j + 1) * 5 + i][k])
            CNOT | (bc[(i + 1) % 5][(k - 1) % 64], state[(j + 2) * 5 + i][k])
            CNOT | (bc[(i + 1) % 5][(k - 1) % 64], state[(j + 3) * 5 + i][k])
            CNOT | (bc[(i + 1) % 5][(k - 1) % 64], state[(j + 4) * 5 + i][k])

    for i in range(5):
        for k in range(64):
            CNOT | (bc[(i + 1) % 5][(k - 1) % 64], state[j * 5 + i][k])
            CNOT | (bc[(i + 4) % 5][k], state[(j + 1) * 5 + i][k])
            CNOT | (bc[(i + 4) % 5][k], state[(j + 2) * 5 + i][k])
            CNOT | (bc[(i + 4) % 5][k], state[(j + 3) * 5 + i][k])
            CNOT | (bc[(i + 4) % 5][k], state[(j + 4) * 5 + i][k])

    for i in range(5):
        for j in range(64):
            CNOT | (state[i][j], bc[i][j])
            CNOT | (state[i + 5][j], bc[i][j])
            CNOT | (state[i + 10][j], bc[i][j])
            CNOT | (state[i + 15][j], bc[i][j])
            CNOT | (state[i + 20][j], bc[i][j])

def phi(eng, state, keccakf_piln, keccakf_rotc):
    t = state[1]
    for i in range(24):
        j = keccakf_piln[i]
        temp = state[j]
        t = ROL(eng, t, keccakf_rotc[i])
        state[j] = t
        t = temp

def chi(eng, copy_state, state, result_state):
    for i in range(5):
        for j in range(5):
            for k in range(64):
                AND_gate_dag(eng, copy_state[i * 5 + ((j + 1) % 5)][k], state[i * 5 + ((j + 2) % 5)][k], result_state[j + i * 5][k])

def chi_reverse(eng, copy_state, state, result_state, anc):
    for i in range(5):
        for j in range(5):
            for k in range(64):
                AND_gate(eng, copy_state[i * 5 + ((j + 1) % 5)][k], state[i * 5 + ((j + 2) % 5)][k], result_state[j + i * 5][k], anc[i][j][k])

print('Estimate cost...')
Resource = ResourceCounter()
eng = MainEngine(Resource)
resource_check = 1
SHA(eng)
print(Resource)
print('\n')
eng.flush()
