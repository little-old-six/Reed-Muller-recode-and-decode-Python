import sys

sys.path.append(r'..\Reed-Muller decoding\Reed-Muller')

from RM_main import RM_main


''' If you can not run this program successfully, try the 'Start.py' in the file of 'Reed-Muller' '''
if __name__ == '__main__':
    R = 2
    M = 5
    N = 3    # number of error

    RM_main(R, M, N)
