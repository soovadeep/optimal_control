function residual = OLcontrolBC(ya,yb)

global x10

residual = [(ya(1) - x10); ya(2); ya(3); ya(4); (yb(5)); (yb(6)); (yb(7)); (yb(8))];