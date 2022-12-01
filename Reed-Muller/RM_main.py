import sys
from RM_recode import RM
from RM_decode import RM_low_radix, RM_high_radix
from Tools import H_compare


def RM_main(R, M, N):
    if R > M:
        sys.exit("Error! R shouldn't more than M")

    # R-M encode
    rm = RM(R, M, N)
    U_0, C = rm.start()

    # R-M decode
    if R == 1:
        a = RM_low_radix(R, M, C)
    else:
        a = RM_high_radix(R, M, C)
    U_mes = a.start()
    print('Decoded message symbol:', U_mes)

    # calculate the bit error rate
    Pe = H_compare(U_0, U_mes) / len(U_mes)
    print("Pe = %.2f" % Pe)
    print('Over!')
