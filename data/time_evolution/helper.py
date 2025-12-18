import numpy as np
## Helper functions


def log(string):
    """
    Print progress
    """

    print("-"*20 + " " + string + " ..." + "-"*20)




def init_state_smart(length, N, sign, path_out, random):
    """
    Function to create initial state txt file:
    """
    print("-"*20 + " Creating initial state (Neel background + holes aranged in stripe order)")

    assert N % 2 ==1

    ## create file
    file_ = open(path_out, 'w')

    afm_bg = [0 for it in range(length)]
    afm_bg[length//2] = 1
    if N<= length: 
        for n in range(N):
            if n<(N-1)/2:
                idx = length//2-n-1
                afm_bg[idx] = (-1)**(idx-length//2)
            elif n>(N-1)/2:
                idx = length//2+n-(N-1)//2
                afm_bg[idx] = (-1)**(idx-length//2)
    print(afm_bg)
    #sign
    afm_bg *= sign

    #write to file
    typ = list(range(0,1))*length

    for i in range(len(afm_bg)):
        file_.write('- *'+ '\n')
        print('- *'+ '\n')
        string = ''
        con = ['0,']
        if afm_bg[i]!=0:
            con[0] = '1,'

            for j in range(1):
                string += con[j]
            if afm_bg[i]>0:
                s="+"
            else:
                s=""
            file_.write('. +1@' +s+ str(0.5*afm_bg[i]) + '\n')
        else:
            for j in range(1):
                string += con[j]

            file_.write('. 1@'+ '0' + '\n')
            print('. 1@0 \n')
