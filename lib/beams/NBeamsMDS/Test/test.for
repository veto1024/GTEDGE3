	program Test

	implicit none
	include 'numbers.inc'

		
	real who(num1), what(num2), why(num3), array(num3,num2)
	real data1(num1), data2(num2), data3(num3)

	namelist/DAT/data1,data2,data3,array

	open(unit=10,file='inputs.dat',status='old')
	read(10,DAT)
	close(unit=10,status='keep')

	who = data1
	what = data2
	why = data3
	
c	array(1,1) = 2.0, 4.0, 6.0, 8.0
c	array(1,2) = 1.0, 3.0, 5.0, 7.0

	open(20,file='ouputs.dat',status='unknown')
	write(20,*) who
	write(20,*) what
	write(20,*) why
	write(20,*) array(1,1:2)
	write(20,*) array(2,1:2)
	write(20,*) array(3,1:2)
	
c	write(20,2400) array(1,1)
c	write(20,2500) array(2,1)
c	write(20,2600) array
	close(20,status='unknown')
	
	end