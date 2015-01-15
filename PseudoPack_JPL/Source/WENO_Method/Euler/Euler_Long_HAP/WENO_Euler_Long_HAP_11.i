        Sum_a = a(0) + a(1) + a(2) + a(3) + a(4) + a(5)

        Omega(0) = a(0)/Sum_a
        Omega(1) = a(1)/Sum_a
        Omega(2) = a(2)/Sum_a
        Omega(3) = a(3)/Sum_a
        Omega(4) = a(4)/Sum_a
        Omega(5) = a(5)/Sum_a

        a(0) = Omega(0)*(d(0)*(d(0)+ONE)+Omega(0)*(Omega(0)-THREE*d(0)))/   &
                        (d(0)* d(0)     +Omega(0)*(ONE     -  TWO*d(0)))
        a(1) = Omega(1)*(d(1)*(d(1)+ONE)+Omega(1)*(Omega(1)-THREE*d(1)))/   &
                        (d(1)* d(1)     +Omega(1)*(ONE     -  TWO*d(1)))
        a(2) = Omega(2)*(d(2)*(d(2)+ONE)+Omega(2)*(Omega(2)-THREE*d(2)))/   &
                        (d(2)* d(2)     +Omega(2)*(ONE     -  TWO*d(2)))
        a(3) = Omega(3)*(d(3)*(d(3)+ONE)+Omega(3)*(Omega(3)-THREE*d(3)))/   &
                        (d(3)* d(3)     +Omega(3)*(ONE     -  TWO*d(3)))
        a(4) = Omega(4)*(d(4)*(d(4)+ONE)+Omega(4)*(Omega(4)-THREE*d(4)))/   &
                        (d(4)* d(4)     +Omega(4)*(ONE     -  TWO*d(4)))
        a(5) = Omega(5)*(d(5)*(d(5)+ONE)+Omega(5)*(Omega(5)-THREE*d(5)))/   &
                        (d(5)* d(5)     +Omega(5)*(ONE     -  TWO*d(5)))
