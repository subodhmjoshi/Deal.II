# Program to solve Ax=b iteratively using method of steepest descent
# Reference: An introduction to the conjugate gradient method without agonizing pain, version 1.25, Jonathan Shewchuck
# This is a demo code thus catering only for the 2d spaces
# Author: Subodh Joshi, IMB (U-Bordeaux) and INRIA BSO

import numpy as np
import matplotlib.pyplot as plt


# matrix A
A = np.array([[3,2],[2,6]]);
# right hand side b
b = np.array([2,-8]);
# constant c for the quadratic form
c = 0.0;


def quadratic_form(x1,x2):
    # x1 and x2 are the x and y vectors created in meshgrid operation.
    # i.e. define your x and y, then create x1,x2 = np.meshgrid(x,y) and pass that as the argument

    # F is the required quadratic form 
    F = np.zeros(shape=(np.shape(x1)[0], np.shape(x2)[1]));
    for i in range(np.shape(x1)[0]):
        for j in range(np.shape(x2)[1]):
            # x1 stores x values and x2 stores y values.
            # extract values to create a single vector x
            vec = np.array([x1[i][j],x2[i][j]]);
            # Quadratic form is given by F = 1/2 x_transpose A x - b_transpose x + c
            # which is same as:
            F[i][j] = 1./2.*np.dot(vec,np.dot(A,vec)) - np.dot(vec,b) + 0.0;
            pass
        pass;
    return F;


def compute_error(x1,x2):
    # Gives L2 norm (distance) between x1 and x2
    return np.sqrt((x1[0] - x2[0])**2 + (x1[1] - x2[1])**2);


def steepest_descent(start):
    # Input:
    # start is the starting point. It should be a 2d vector storing x and y coordinates of the starting point.

    # Output:
    # P, X, Y, E, Iteration described as follows:
    
    # P stores the points on (x,y) plane. Note: P[i] = (X[i], Y[i])
    P = [];
    P.append(start);

    # X and Y arrays store the x and y coordinates of the subsequent steps 
    X = [];
    X.append(start[0]);
    Y = [];
    Y.append(start[1]);

    # E stores the global error 
    E = [];

    error = 10**10; # arbitrary large number
    Tol = 10**-3;   # small number

    Iteration = 0;

    while error > Tol:
        # REMEMBER residual is the direction of the steepest descent!
        residual = b - np.dot(A,P[-1]);
        
        # alpha is the distance along r for the next step 
        # alpha is given by (r_transpose * r)/(r_transpose * A * r) which is same as:
        alpha = np.dot(residual,residual)/ float(np.dot(residual,np.dot(A,residual)));
        
        # update:
        next_step = P[-1] + alpha * residual;

        error = compute_error(next_step, P[-1]);

        # storing data
        E.append(error);
        P.append(next_step);
        X.append(P[-1][0]);
        Y.append(P[-1][1]);
        Iteration += 1;
        pass;

    return P, X, Y, E, Iteration;


def plotC(X,Y,F):
    # plots the contour plots of the quadratic form
    CS = plt.contour(X,Y,F,10);
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xlabel('x1')
    plt.ylabel('x2')
    plt.grid();



if __name__=="__main__":
    # Create data given in the document
    x1 = np.linspace(-4,6,50);
    x2 = np.linspace(-6,4,50);
    X,Y = np.meshgrid(x1,x2)

    F = quadratic_form(X,Y);
    plt.figure
    plotC(X,Y,F);

    # Initial guess given in the document
    start = np.array([-2,-2])
    P,X1,Y1,E,Iteration = steepest_descent(start); 

    plt.plot(X1,Y1);
    plt.title("Total number of iterations: "+str(Iteration));
    plt.show();

    plt.figure
    plt.semilogy(np.arange(Iteration),E);
    plt.xlabel("Iterations")
    plt.ylabel("Error");
    plt.show()




