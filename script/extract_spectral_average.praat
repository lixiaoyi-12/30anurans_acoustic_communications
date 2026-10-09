saveFileName$="~/Downloads/postdoc工作/1海南蛙类鸣声与ABR/鸣声数据/all/Rhacophorus_rhodopus红蹼树蛙/红蹼树蛙7_spectral_average.txt"

binNum=Get number of bins
binW=Get bin width

fileappend 'saveFileName$' Hz
fileappend 'saveFileName$' ,
fileappend 'saveFileName$' dB
fileappend 'saveFileName$' 'newline$'

for n from 1 to 40

    targetF=200*n
    lowF=targetF-100
    highF=targetF+100

    sumPower=0
    count=0

    for j from 1 to binNum

        pitchV=(j-1)*binW

        if pitchV>=lowF
            if pitchV<highF

                rV=Get real value in bin... 'j'
                iV=Get imaginary value in bin... 'j'

                power=rV*rV+iV*iV

                sumPower=sumPower+power
                count=count+1

            endif
        endif

    endfor

    if count>0

        meanPower=sumPower/count

        dB=10*log10(meanPower/(4*10^-10))+15

        fileappend 'saveFileName$' 'targetF'
        fileappend 'saveFileName$' ,
        fileappend 'saveFileName$' 'dB'
        fileappend 'saveFileName$' 'newline$'

    endif

endfor