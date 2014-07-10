make

if [ $? -ne 0 ]
then
	exit 1
fi

#8000 so that the arrays can fit inside the texture memory ( 8 * 8000 = 64000 < 65000 )
for j in 1000 8000
do

	echo "$j"

	for i in sog
	do
		./smatcher $i -m 8 -p_size $j -n 3999744 -alphabet 2
	done

	echo ""
	
	for i in sog
	do
		./smatcher $i -m 8 -p_size $j -n 4628736 -alphabet 4
	done

	echo ""
	
	for i in sog
	do
		./smatcher $i -m 8 -p_size $j -n 116234496 -alphabet 4
	done

	echo ""	
	
	for i in sog
	do
		./smatcher $i -m 8 -p_size $j -n 177649920 -alphabet 20
	done

	echo ""
	
	for i in sog
	do
		./smatcher $i -m 8 -p_size $j -n 10821888 -alphabet 20
	done

	echo ""

	for i in sog
	do
		./smatcher $i -m 8 -p_size $j -n 1903104 -alphabet 128
	done

	echo ""

done
