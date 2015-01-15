        Sum_a = a(0) + a(1) + a(2) + a(3) 

        Omega(0) = a(0)/Sum_a
        Omega(1) = a(1)/Sum_a
        Omega(2) = a(2)/Sum_a
        Omega(3) = a(3)/Sum_a

        a(0) = Omega(0)*(d(0)*(d(0)+ONE)+Omega(0)*(Omega(0)-THREE*d(0)))/   &
                        (d(0)* d(0)     +Omega(0)*(ONE     -  TWO*d(0)))
        a(1) = Omega(1)*(d(1)*(d(1)+ONE)+Omega(1)*(Omega(1)-THREE*d(1)))/   &
                        (d(1)* d(1)     +Omega(1)*(ONE     -  TWO*d(1)))
        a(2) = Omega(2)*(d(2)*(d(2)+ONE)+Omega(2)*(Omega(2)-THREE*d(2)))/   &
                        (d(2)* d(2)     +Omega(2)*(ONE     -  TWO*d(2)))
        a(3) = Omega(3)*(d(3)*(d(3)+ONE)+Omega(3)*(Omega(3)-THREE*d(3)))/   &
                        (d(3)* d(3)     +Omega(3)*(ONE     -  TWO*d(3)))
