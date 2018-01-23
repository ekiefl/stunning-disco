import numpy as np
def linear_regression(Xm,y,theta,T,alpha):
    """
    Performs multi-variate linear regression.
    
    Fits n-variate data x to
    
    h(x;theta) = theta0 + theta1*x1 + ... + thetan*xn
    
    by using method of gradient descent to find the theta 
    that minimizes the least mean squared error cost function
    
    J(theta) = 1/(2m) * sum_i( [h(x^(i);theta) - y^(i)]^2 )
    
    where m is the number of data points and x^(i) is the ith
    data points.
    
    (https://en.wikipedia.org/wiki/Linear_regression)
    
    INPUTS
    ------
    Xm : numpy array
        Xm is a matrix of explanatory variables with 
        shape = (a,b), where a is the number of observations
        and b is the number of explanatory variables. In the
        case of a 2D linear regression, the matrix should
        look like either
        
        var1 var2
        [ 9   9 ]
        [ 8   2 ]
        [ 3   2 ]
        [... ...]
        
        or
        
        var0 var1 var2
        [ 1   9   9 ]
        [ 1   8   2 ]
        [ 1   3   2 ]
        [... ... ...]
        
        In the second case the first row is the 
        pseudo explanatory variable for the parameter theta0
        and is necessarily all 1's. If the user provides
        Xm like the first type, this column will be
        automatically added. If you pass Xm with the column
        already added, it has to be column 0.
    y : numpy array
        y is the dependent variable and should have a length
        equal to the number of Xm rows.
    theta : numpy array
        y is the initial guess for your parameter estimates.
        If you're data is n-dimensional, theta should be of 
        length n+1 and ordered. For example, if n=2, 
        theta = [theta0,theta1,theta2]. I'm being
        redundant at this point, but just to drive home the 
        point, the fit function would then be
        y = theta0 + theta1*x1 + theta2*x2.
    T : int
        The maximum number of iterations performed. If 
        convergence is not reached by T iterations, the params
        are returned.
    alpha : float
        The step size of the gradient descent. alpha should be
        large enough so that theta traverses the parameter 
        landscape sufficiently, but small enough so that the
        cost function J decreases with each iteration.
    conv : float, default = None
        This value defines convergence of the cost function.
        throughout each iteration t, if cost[t+1]-cost[t] < conv,
        the fit is said to have converged. It is situation
        specific and if not provided, all T iterations
        will be carried out.
    
    RETURNS
    -------
    theta_history : numpy array
        numpy array of shape (iters,n+1), where iters
        is the number of iterations and n is the dimension of the
        data. Each theta_history[t,:] defines theta for iteration t
    cost_history : numpy array
        1D numpy array of length iters, where iters
        is the number of iterations. Each cost_history[t] defines
        the cost at each iteration t.
    """
    
#   append pseudo-column if not present
    if (Xm.ndim == 1):
        pseudo = np.ones(np.shape(Xm)[0])
        Xm = np.vstack((pseudo,Xm.T)).T
    elif all(Xm[:,0]!=1):
        pseudo = np.ones(np.shape(Xm)[0])
        Xm = np.vstack((pseudo,Xm.T)).T
    else:
        pass
    
        
#   feature scaling with mean normalization so J is not highly 
#   elliptical and exists about 0 vector the means and range are
#   stored so theta can be converted back once the fit converges
    ranges = np.zeros(np.shape(Xm)[1]-1)
    means = np.zeros(np.shape(Xm)[1]-1)
    rangey = np.max(y)-np.min(y)
    meany = np.mean(y)
    y = (y-meany)/rangey

    for col in range(len(ranges)):
        ranges[col] = np.max(Xm[:,col+1])-np.min(Xm[:,col+1])
        means[col] = np.mean(Xm[:,col+1])
        Xm[:,col+1] = (Xm[:,col+1]-means[col]) / ranges[col]

#   track the cost and parameters
    cost_history = np.zeros(T)
    theta_history = np.zeros((T,2))
    
#   iterate 
    for t in range(0,T):
        cost_history[t] = calc_cost(Xm,y,theta)
        theta_history[t,:] = theta
        theta = take_step(Xm,y,theta,alpha)
        
#       stop if it converges
        if conv != None:
            if (cost_history[t-1]-cost_history[t]<conv) & (t > 0):
                theta_history[t,:] = theta
                cost_history[t] = calc_cost(Xm,y,theta)
                theta_history = theta_history[cost_history!=0,:]
                cost_history = cost_history[cost_history!=0]
                break
    
#   reverse effects of feature scaling and mean normalization
    nasty = np.sum(theta_history[:,1:]*(means/ranges),axis=1)
    theta_history[:,0] = (theta_history[:,0]-nasty)*rangey + meany
    theta_history[:,1:] = rangey*theta_history[:,1:]/ranges
    
    return cost_history, theta_history

def take_step(Xm,Y,theta,alpha):
    
#   hypothesis
    h = np.dot(Xm,theta)
#   gradient
    grad = 1./m * np.dot(Xm.T,(h-Y))
#   each step is just -grad * alpha
    step = -grad * alpha
#   return the new theta
    return theta + step

def calc_cost(Xm,y,theta):
#   hypothesis
    h = np.dot(Xm,theta)
#   sum of squares
    cost = 1./(2*m) * np.sum((h-y)**2)
    return cost

##################################################
