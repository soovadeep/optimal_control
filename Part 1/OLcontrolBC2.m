function residual = OLcontrolBC2(ya,yb)

global x10 

residual = [(ya(1) - x10); ya(2); ya(3); ya(4); (yb(1)); (yb(2)); (yb(3)); (yb(4))];