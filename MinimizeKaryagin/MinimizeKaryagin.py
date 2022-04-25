import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


EPS = 0.000001


#exercise 3
def OriginalFunction(x):
    return np.power(x,2) - 2*x - 2*np.cos(x)

def DerivativeFunction(x):
    return 2*x-2+2*np.sin(x)

def SecondDerivativeFunction(x):
    return 2+2*np.cos(x)

#3.1 
def Plotting(leftBorderX, rightBorderX,x):
    vectorX = np.arange(-0.5, 1, 0.0000001)
    plt.plot(vectorX, OriginalFunction(vectorX), x, OriginalFunction(x))
    plt.show()

def BitwiseSearchMethod(leftBorderX, rightBorderX):
    stepCount = 0;
    step = 4*EPS
    x = leftBorderX
    pointsOfApproach = []
    while np.abs(step) > EPS:
       if OriginalFunction(x) < OriginalFunction(x + step):
            step = -(step/4)
            x += step
            pointsOfApproach.append(x)
       else: 
            x+=step
            pointsOfApproach.append(x)
       stepCount+=1;
    print("Number of steps of bitwise search method:" , stepCount)
    Plotting( leftBorderX, rightBorderX, np.array(pointsOfApproach) )
    return x


#3.2
def DichotomyMethod(leftBorderX, rightBorderX):
    stepCount = 0
    leftBorderForPlot = leftBorderX
    rightBorderForPlot = rightBorderX
    pointsOfApproach = []
    epsN = (rightBorderX - leftBorderX) / 2
    while epsN > EPS:
        firstX = (rightBorderX + leftBorderX - 1.5*EPS)/2
        pointsOfApproach.append(firstX)
        secondX = (rightBorderX + leftBorderX + 1.5*EPS)/2
        pointsOfApproach.append(secondX)
        if OriginalFunction(firstX) < OriginalFunction(secondX):
            rightBorderX = secondX
        else:
            leftBorderX = firstX
        epsN = (rightBorderX - leftBorderX) / 2
        stepCount+=1
    print("Number of steps of dichotomy method:" , stepCount)
    resaultX = (rightBorderX + leftBorderX)/2
    pointsOfApproach.append(resaultX)
    Plotting( leftBorderForPlot, rightBorderForPlot, np.array(pointsOfApproach) )
    return resaultX


#3.3
def FindX2(x1,x3):
    x2 =x1+0.05
    while OriginalFunction(x1) < OriginalFunction(x2) or OriginalFunction(x3) < OriginalFunction(x2):
        x2+=0.05
    return x2

def ParabolaMethod(leftBorderX, rightBorderX):
    x1 = leftBorderX
    x3 = rightBorderX
    x2 = FindX2(x1,x3)
    pastX = 10;
    resaultX = 0
    stepCount = 0;
    pointsOfApproach = []

    while np.abs(pastX-resaultX) > EPS:
        a1 = (OriginalFunction(x2) - OriginalFunction(x1))/(x2-x1)
        a2 = ((OriginalFunction(x3)-OriginalFunction(x1))/(x3-x1) - a1) * 1/(x3-x2)
        pastX = resaultX
        resaultX = 0.5*(x1+x2-a1/a2)
        if resaultX>x2:
            x1 = x2
            x2 = resaultX
        else:
            x1 = resaultX
        stepCount+=1
        pointsOfApproach.append(resaultX)
    print("Number of steps of parabola method:" , stepCount)
    Plotting( leftBorderX, rightBorderX, np.array(pointsOfApproach) )
    return resaultX


#3.4
def NewtonMethod(leftBorderX, rightBorderX):
    x = (rightBorderX + leftBorderX) / 2
    while np.abs(DerivativeFunction(x)) > EPS:
        x = x -  DerivativeFunction(x) / SecondDerivativeFunction(x)
    return x

    
#7
def FunctionForDichot(x):
    return np.exp(x)-1-x-np.power(x,2)/2-np.power(x,3)/3

def DichotomyMethodFor7Part(leftBorderX, rightBorderX):
    stepCount = 0
    leftBorderForPlot = leftBorderX
    rightBorderForPlot = rightBorderX
    EPS = 0.0001
    epsN = (rightBorderX - leftBorderX) / 2

    while epsN > EPS:
        firstX = (rightBorderX + leftBorderX - 1.5*EPS)/2
        secondX = (rightBorderX + leftBorderX + 1.5*EPS)/2
        if FunctionForDichot(firstX) < FunctionForDichot(secondX):
            rightBorderX = secondX
        else:
            leftBorderX = firstX
        epsN = (rightBorderX - leftBorderX) / 2
        stepCount+=1
    print("Number of steps of dichotomy method:" , stepCount)
    resaultX = (rightBorderX + leftBorderX)/2
    return resaultX


#8
def FunctionForNewton(x):
    return x*np.arctan(x)-0.5*np.log(1+np.power(x,2))

def DerivativeFunctionForNewton(x):
    return np.arctan(x)

def SecondDerivativeFunctionForNewton(x):
    return 1/(np.power(x,2)+1)

def NewtonMethodFor8Part():
    rightBorderX = -3
    leftBorderX = 5
    x = (rightBorderX + leftBorderX) / 2
    while np.abs(DerivativeFunctionForNewton(x)) > EPS:
        x = x -  DerivativeFunctionForNewton(x) / SecondDerivativeFunctionForNewton(x)
    return x


def TauCalculating(x):
    return np.power(DerivativeFunction(x),2)/ (np.power( DerivativeFunction(x) , 2) + np.power( DerivativeFunction(x - DerivativeFunction(x)/SecondDerivativeFunction(x)),2))

def NewtonRaphsonMethod(leftBorderX, rightBorderX):
    x = (rightBorderX - leftBorderX) / 2
    while np.abs(DerivativeFunctionForNewton(x)) > EPS:
        x = x - TauCalculating(x) * DerivativeFunctionForNewton(x) / SecondDerivativeFunctionForNewton(x)
    return x


def main():
    print("Golden function resault:", optimize.golden(OriginalFunction, brack=(-0.5,1), full_output=False))
    print("fminbound function resault:", optimize.fminbound( OriginalFunction, -0.5, 1))
   

    print("bitwise search method resault: " , BitwiseSearchMethod(-0.5,1))
    print("of dichotomy method resault:", DichotomyMethod(-0.5,1) )
    #print("Number of steps of midpoint method: " , MidpointMethod(-0.5,1))
    print("parabola method resault: " , ParabolaMethod(-0.5,1))
    print("Newton method resault: " , NewtonMethod(-0.5,1) )

    print("Dichotomy method resault:", DichotomyMethodFor7Part(-5,5) )

    print("Newton method resault:", NewtonMethodFor8Part())
    print("Number of steps of Newton - Raphson method: " , NewtonRaphsonMethod(-0.5,1) )

    return 0

main()