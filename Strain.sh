#!/bin/bash



# This scripts is built for finding optimal relaxed point, when strain is applied.

export PATH=/THL4/home/haiyan/ziang/bin:$PATH

# VASP bin

BIN="yhrun -pIOR -N2 -n24 /THL4/home/haiyan/bin/vasp_ncl"

# BIN="mpirun -n 1 vasp_ncl"

# BIN="/home/allegro/vasp-openblas/bin/mpirun -np 12 /home/allegro/vasp-openblas/vasp.6.1.0/bin/vasp_ncl"

# BIN="/opt/nvidia/hpc_sdk/Linux_x86_64/20.9/comm_libs/mpi/bin/mpirun -n 1 /home/allegro/vasp-openacc/vasp.6.1.0/bin/vasp_ncl"



# set strain and it's applied function

mode=0 #set 0 for biaxial strain, set 1 for uniaxial strain

strain=0.03 # strain=da/a = db/b for mode0, strain=dc/c for mode1



# set the iteration parameter

itenumber=10 

initpercent=-0.03

teststep=0.1

# method of bisection entry threshold

# sigma=0.05



# set the lattice constant

a0=`nl POSCAR | sed -n '3p' | awk '{print $2}'`

b0=`nl POSCAR | sed -n '4p' | awk '{print $3}'`

c0=`nl POSCAR | sed -n '5p' | awk '{print $4}'`



# Calculate module, It's output is the energy without entropy.

Runcalculate()

{

    # enter calculate direction

    mkdir test$test

    cp -t test$test INCAR POTCAR KPOINTS POSCAR

    cd test$test

    # build POSCAR

    sed -i "3c $1      0      0" POSCAR

    sed -i "4c 0      $2      0" POSCAR

    sed -i "5c 0      0      $3" POSCAR



    $BIN >> ../runtime

    if grep "reached required accuracy" OUTCAR >> ../runtime

    then

        echo "`date` Test$test is completed" >> ../log

        echo "$test  $1  $2  $3  `cat OUTCAR | grep "energy without entropy" | tail -n1`  `tail -n1 OSZICAR`" >> ../data

        echo `cat OUTCAR | grep "energy without entropy" | tail -n1 | awk '{print $5}'`

    else

        echo "`date` Error has been appeared, when processing Test$test" >> ../log

        echo 1

    fi

    cd ..    

}



#initial output file

echo "Index | a            | b                  | c                   | Energy without entropy                                                             | OSZICAR " >> data



#initial POSCAR parameter

a=`echo "(1+(1-$mode)*$strain+$mode*$initpercent) *$a0" | bc`

#a=!b, try using radio to get b in uniaxial mode.

#b=`echo "(1+(1-$mode)*$strain+$mode*$initpercent) *$b0" | bc`

b=`echo "(1+(1-$mode)*$strain+$mode*($a/$a0-1)) *$b0" | bc`

c=`echo "(1+$mode*$strain+(1-$mode)*$initpercent) *$c0" | bc`



#initial optimize algorithm, you can use external parameter to replace it.

test=0

ta=($a 0 $a 0 $a)

tb=($b 0 $b 0 $b)

tc=($c 0 $c 0 $c)

e=(0 0 0 0 0)

while getopts "t:a:b:c:e:" arg

do

    case $arg in

        t)

          test=$OPTARG

          ;;

        a)

          index=0

          for argu in $OPTARG

          do

            ta[$index]=$argu

            index=$[$index+1]

          done

          ;;

        b)

          index=0

          for argu in $OPTARG

          do

            tb[$index]=$argu

            index=$[$index+1]

          done

          ;;

        c)

          index=0

          for argu in $OPTARG

          do

            tc[$index]=$argu

            index=$[$index+1]

          done

          ;;

        e)

          index=0

          for argu in $OPTARG

          do

            e[$index]=$argu

            index=$[$index+1]

          done

          ;;

        ?)

          echo -e "-t test number; -a ta array; -b tb array; -c tc array; -e e array\nNote: You should  use string to input five elements array. like this, \"1 2 3 4 5\""

          exit 1

          ;;

    esac

done



echo "test $test | ta ${ta[*]} | tb ${tb[*]} | tc ${tc[*]} | e ${e[*]}"

k1=1

k2=1

iteration=0

while [ $iteration -lt $itenumber ]

do

#build *1,*2,*3 temp parameter for slope cauculate, depend on teststep that set above.

if [ $test -lt 3 ]

then

    ta[$[$test * 2]]=`echo "$mode*($test-1)*$teststep+${ta[$[$test * 2]]}" | bc`

    tb[$[$test * 2]]=`echo "$mode*(${ta[$[$test * 2]]}/$a-1)*${tb[$[$test * 2]]}+${tb[$[$test * 2]]}" | bc`

    tc[$[$test * 2]]=`echo "(1-$mode)*($test-1)*$teststep+${tc[$[$test * 2]]}" | bc`

    e[$[$test * 2]]=`Runcalculate ${ta[$[$test * 2]]} ${tb[$[$test * 2]]} ${tc[$[$test * 2]]}`

    test=$[$test+1]

else

    if [ $mode -eq 0 ]

    then

        k1=`echo "scale=8; (${e[2]} - ${e[0]}) / (${tc[2]} - ${tc[0]})" | bc`

        k2=`echo "scale=8; (${e[4]} - ${e[2]}) / (${tc[4]} - ${tc[2]})" | bc`

        tag=`echo "scale=8; $k1 * $k2" | bc`

        # if [ $tag -lt 0 ]

        if [ $(echo "$tag < 0" | bc) -eq 1 ] 

        then                       

            ta[1]=${ta[0]}

            tb[1]=${tb[0]}

            tc[1]=`echo "scale=8; (${tc[0]} + ${tc[2]})/2" | bc`

            ta[3]=${ta[0]}

            tb[3]=${tb[0]}

            tc[3]=`echo "scale=8; (${tc[2]} + ${tc[4]})/2" | bc`

            e[1]=`Runcalculate ${ta[1]} ${tb[1]} ${tc[1]}`

            test=$[$test+1]

            e[3]=`Runcalculate ${ta[3]} ${tb[3]} ${tc[3]}`

            test=$[$test+1]

            run=0

            tag=1

            # while [ $tag -gt 0 ]

            while [ $(echo "$tag > 0" | bc) -eq 1 -a $run -lt 3 ]

            do

                k1=`echo "scale=8; (${e[$[$run + 1]]} - ${e[$[$run + 0]]}) / (${tc[$[$run + 1]]} - ${tc[$[$run + 0]]})" | bc`

                k2=`echo "scale=8; (${e[$[$run + 2]]} - ${e[$[$run + 1]]}) / (${tc[$[$run + 2]]} - ${tc[$[$run + 1]]})" | bc`

                tag=`echo "scale=8; $k1 * $k2" | bc`

                run=$[$run+1]

                # if the relaxation approach to maximum value, switching to next circulation.

                if [ $(echo "$k1 > 0" | bc) -eq 1 ]

                then

                    tag=1

                fi

            done

            # use 402 assignment to avoiding override

            tc[4]=${tc[$[$run + 1]]}

            tc[0]=${tc[$[$run - 1]]}

            tc[2]=${tc[$[$run + 0]]}            

            e[4]=${e[$[$run + 1]]}

            e[0]=${e[$[$run - 1]]}

            e[2]=${e[$[$run + 0]]}

            # if the maximum appeared at run=3, expland test section for optimization.

            if [ $tag -eq 1 -a $run -eq 3 ]

            then

                # Left shift

                tc[4]=${tc[2]}

                tc[2]=${tc[0]}

                tc[0]=`echo "(-1)*(1-$mode)*$teststep+${tc[0]}" | bc`

                e[4]=${e[2]}

                e[2]=${e[0]}

                e[0]=`Runcalculate ${ta[0]} ${tb[0]} ${tc[0]}`

                test=$[$test+1]

            fi 

        else

            # if [ $k1 -lt 0 ]

            if [ $(echo "$k1 < 0" | bc) -eq 1 ]

            then

                # Right shift

                tc[0]=${tc[2]}

                tc[2]=${tc[4]}

                tc[4]=`echo "(1-$mode)*$teststep+${tc[4]}" | bc`

                e[0]=${e[2]}

                e[2]=${e[4]}

                e[4]=`Runcalculate ${ta[4]} ${tb[4]} ${tc[4]}`

                test=$[$test+1]

            else

                # Left shift

                tc[4]=${tc[2]}

                tc[2]=${tc[0]}

                tc[0]=`echo "(-1)*(1-$mode)*$teststep+${tc[0]}" | bc`

                e[4]=${e[2]}

                e[2]=${e[0]}

                e[0]=`Runcalculate ${ta[0]} ${tb[0]} ${tc[0]}`

                test=$[$test+1]

            fi

        fi

    else

        echo "Uniaxial strain calculate have not being present yet."

    fi

fi

# e1=`Runcalculate $ta1 $tb1 $tc1`

echo "Test $[$test - 1] | ta ${ta[*]} | tb ${tb[*]} | tc ${tc[*]} | e ${e[*]} | k1 $k1 | k2 $k2 | run $run" >> log

echo -e "-t $test\n-a ${ta[*]}\n-b ${tb[*]}\n-c ${tc[*]}\n-e ${e[*]}" > resume

iteration=$[$iteration+1]

done
