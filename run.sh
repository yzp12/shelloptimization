a=1.5
b=1.7
c=`echo "$a + $b"|bc`
echo $c
for((i=1;i<=10;i++));do shelloptimization `echo "$i * 0.01"|bc`;done

