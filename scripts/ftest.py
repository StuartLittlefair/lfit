#!/scisoft/i386/bin/python
from __future__ import division
from __future__ import print_function
from builtins import input
from past.utils import old_div
import scipy.special
import sys

if len(sys.argv) == 5:
    chi1 = float(sys.argv[1]) 
    dof1 = float(sys.argv[2]) 
    chi2 = float(sys.argv[3]) 
    dof2 = float(sys.argv[4]) 
else:
    chi1 = float(input('> Give chisq for simple model: '))
    dof1 = float(input('> Give D.O.F for simple model: '))
    chi2 = float(input('> Give chisq for complex model: '))
    dof2 = float(input('> Give D.O.F for simple model: '))

val = (chi1-chi2)/(dof1-dof2)/(old_div(chi2,dof2))
ftest = scipy.special.betainc(old_div(dof2,2.),old_div((dof1-dof2),2.),(old_div(dof2,(dof2 + (dof1-dof2)*val))))

print("\n\nF-test value is ", val)
print("Larger values imply better fits with the more complex model")
print("The probability that a random dataset might give a F-test value")
print("exceeding this when the data is actually well represented by the simple model")
print("is %7.4f" % (100.0*ftest), '%')
print("")
print("So you can reject the null hypothesis that the complex model is no better than the simple")
print("model at the %7.4f" % (100.0*(1.0-ftest)), "% confidence level")
