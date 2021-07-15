
with open('FORCE_CONSTANTS','r') as fin:
    with open('out','w') as fout:

        for ii in range(2):
            for jj in range(16):
                
                inds = fin.readline().strip().split()
                inds = [int(x)-1 for x in inds]           
                if inds[0] <= 7: 
                    tmp_1 = 0
                else:
                    tmp_1 = 1
                if inds[1] <= 7:
                    tmp_2 = 0
                else:
                    tmp_2 = 1
                fout.write(f' {tmp_1} {tmp_2} {inds[0]} {inds[1]}\n') # {inds[0]} {inds[1]}\n')

                for kk in range(3):
                    fc = fin.readline().strip().split()
                    fc = [float(x) for x in fc]
                    fout.write(f'{fc[0]: 2.6f} {fc[1]: 2.6f} {fc[2]: 2.6f}\n')




