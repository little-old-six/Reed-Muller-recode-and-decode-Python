

def H_compare(M0, M1):
    dis = 0
    l = len(M0)
    for i in range(l):
        if M1[i] != M0[i]:
            dis += 1
    return dis
