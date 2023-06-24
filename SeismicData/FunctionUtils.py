


def abs_maxF(f,x0,x1,n_points):
    # f is a lambda scalar valued function defined in [x0,x1]

    step   = (x1-x0)/n_points
    x_eval = (xi*step for xi in range(0,n_points+1)) #Evaluation test points
    #x_eval = [xi*step for xi in range(0,n_points+1)]
    return max([abs(f(x)) for x in x_eval])
    
def fplot(ax,f,x0,x1,translation=0,n=5):
    # plots a lambda real function f in a defined interval [x0,x1] in n line segments, in a given axis object.
    
    step = (x1-x0)/n
    x_eval = [x0+xi*step for xi in range(0,n+1)] #Evaluation test points
    f_eval = [f(x) for x in x_eval]

    x_plot = [x+translation for x in x_eval]

    ax.plot(x_plot,f_eval,'b') 
    return ax

def plot_picewise_lambda(ax,F,X,N=5):
    n = len(F) # Extract interval number (lenX)=n+1)

    for t in range(n):
        ax = fplot(ax,F[t],X[t],X[t+1],0, N)

    return ax
     
def plot_picewise_lambda_localX(ax,F,X, N=5):
    n = len(F) # Extract interval number (lenX)=n+1)
    
    for t in range(n):
        ax = fplot(ax,F[t],0,X[t+1]-X[t],X[t],N)

    return ax