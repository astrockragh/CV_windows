import numpy as np
import matplotlib.pyplot as plt

def hist_fit(function, data, fit_type = 'chi2', bound = None, bins = 40, ax = None, print_level = 0, guesses = None):
    """
    This function takes a function or a string of a function. Currently supported: Gauss, Uniform and Poisson. 
    The function take following parameters:
    ---------------------------------------------
    function: The function of the fit
    data:     Data to put in hist and fit
    bound:    What the range of the histogram should be. If none given the min and max of the data is used
    fit_type: How to fit. Availible options are "Chi2", "bllh" and "ullh"
    bins:     Amount of binds. This is ignored if fit_type == 'ullh'
    ax:       What axes the histogram should be painted on.

    ----------------------------------------------
    resturns: 
    migrad_object form the fit
    """
    from scipy.stats import poisson, norm

    # Definition of some typical fits that I would like to call with string.
    fit_functions = {'gauss': lambda x, A, mu, sigma: A  *norm.pdf(x, mu, sigma),\
                 'uniform': lambda x, A: A ,\
                 'poisson': lambda x, A, Lamb: A * norm.pmf(x, Lamb)}

    fit_guesses = {'gauss': lambda data: {'A': len(data), 'mu': np.mean(data), 'sigma': np.std(data)}, \
               'uniform': lambda data: {'A': np.sum(data)/len(data)}}

    fit_types = ['chi2', 'bllh', 'ullh']

    from ExternalFunctions import Chi2Regression, BinnedLH, UnbinnedLH, add_text_to_ax, nice_string_output
    from scipy.stats import chi2
    from iminuit import Minuit
    if type(function) == str:
            
        if function in fit_guesses.keys():
            guesses = fit_guesses[function](data)
        try:
            function = fit_functions[function]
        except KeyError:
            raise NotImplementedError
        

    if fit_type not in fit_types:
        raise NotImplementedError

    if fit_type in ['chi2', 'bllh']:
        if not bound:
            bound = (min(data), max(data))
        y, bin_edges = np.histogram(data, bins = bins, range = bound)
        x = (bin_edges[1:] + bin_edges[:-1])/2
        binwidth = (bound[1] - bound[0])/bins
        sy = np.sqrt(y)
        if fit_type == 'chi2':
            x, y, sy = x[y > 0], y[y > 0], sy[y > 0]
            fit_object = Chi2Regression(function, x, y / binwidth, sy = sy / binwidth)
        else:
            fit_object = BinnedLH(function, data, bins = bins, bound = bound, extended = True)
    else:
        if ax:
            if not bound:
                bound = (min(data), max(data))
            y, bin_edges = np.histogram(data, bins = bins, range = bound)
            x = (bin_edges[1:] + bin_edges[:-1])/2
            binwidth = (bound[1] - bound[0])/bins

        fit_object = UnbinnedLH(function, data, extended = True)
    
    if print_level != 1:
        print_level = 0
    

    if guesses:
        guesses = ", ".join([f"{i} = {j}" for i, j in zip(guesses.keys(), guesses.values())])
        fit_min = eval(f"Minuit(fit_object, pedantic = False, print_level = print_level, {guesses})")
    else:
        fit_min = Minuit(fit_object, pedantic = False, print_level = print_level)
    
    fit_min.migrad()


    if ax:
        ax.plot(x, y, 'r.')
        ax.errorbar(x, y, np.sqrt(y), ls = 'none', elinewidth = 1, capsize = 2, c= 'k')
        x_fit = np.linspace(bound[0], bound[1], 1000)
        y_fit = function(x_fit, *fit_min.args)
        ax.plot(x_fit, binwidth * y_fit, 'k')
        ax.set_ylabel(f"freq/{binwidth:.3f}")

    return fit_min

def weighted_mean(x, errs, ax = None, coords = (0.1, 0.9), dec = 3):
    """
    This function takes as input measurents and errors and returns the weighted mean along with the error.
    The weighted mean is calculated by doing a Chi-Square fit with a constant.
    if ax is given, the function will plot the fit, data and errors on it. 

    The function return weighted_mean, err, and a dictionairy with chi-square value and p-value
    """
    def constant(x, k): return k
    from ExternalFunctions import Chi2Regression, nice_string_output, add_text_to_ax
    from iminuit import Minuit
    from scipy.stats import chi2
    
    ticks = np.arange(len(x))

    chi2_object = Chi2Regression(constant, ticks, x, errs)
    chi2_min = Minuit(chi2_object, pedantic = False)
    chi2_min.migrad()

    Chi2 = chi2_min.fval
    p_val = chi2.sf(Chi2, len(x)-1)
    k = chi2_min.args[0]
    err = chi2_min.errors[0]


    if ax:
        ax.plot(ticks, x, 'r.')
        ax.errorbar(ticks, x, errs, c = 'k', elinewidth = 1, \
                       capsize = 2, ls = 'none')
        ax.hlines(k, min(ticks), max(ticks), ls = '--', label = "Weighted mean")
        d = {"chisquare":   Chi2, \
             "p:":           p_val, \
             "mean:":        k,\
             "mean_err":     err}

        add_text_to_ax(*coords, nice_string_output(d, decimals=dec), ax, fontsize = 10)
        ax.legend()

    return k, err, {"chi2": Chi2, "p": p_val}

def data_display(data, chauvenets = False, significance = 0.5, verbose = False, fig = None, ax = None):
    """
    Display of data.
    Possible use of chauvenets criteria to remove outliers with given significance
    """
    
    # Imports 
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.stats import norm

    # Setup the plots
    if not fig and not ax:
        fig, ax = plt.subplots(figsize = (14, 8))
    
    divider = make_axes_locatable(ax)
    ax1 = divider.append_axes("left", size = "45%", pad = 0.40) # Divide plots
    
    # Calculate mean and std
    data = np.asarray(data)
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    
    # Plot params
    xmin, xmax = min(data) -std, max(data) + std
    xs = np.linspace(xmin, xmax, 1000)
    
    # Setup a mask
    mask = np.ones(len(data)).astype(bool)
    
    if chauvenets:
        # TO DO: Change plot in ax1 to scatter. To show rejected points.
        
        # First plot with the primary measurements before applying
        ax1.errorbar(np.arange(len(data)), data, std, elinewidth = 1, capsize = 2, ls = 'none', c = 'k', alpha = 0.1)
    
        ax1.hlines(mean, -0.5, len(data) + 0.5, alpha = 0.2)
        ax1.hlines((mean - std, mean + std), -0.5, len(data) + 0.5, ls = '--', alpha = 0.2)
        
        ax.plot(xs, norm.pdf(xs, mean, std), 'k--', alpha = 0.35)
        ax.scatter(data, norm.pdf(data, mean, std), c = 'gray', zorder = 5, label = "Before Chauvenets", alpha = 0.5)

        # Copy the data
        data_copy = data.copy()
        mean_copy, std_copy = mean, std
        N=len(data)
        k=0
        while True:
            # Loop: Calculate pull -> p - value
            zs = abs(data_copy - mean_copy)/std_copy
            zs[np.logical_not(mask)] = 0 
            index = np.argmax(zs)
            p = 2 * norm.sf(zs[index]) 
            accept = 1-(1-p)**(N-k) > significance
            k+=1
            if verbose:
                print(f"mean: {mean_copy:.4f}, \t RMS: {std_copy:.4f}, \t point: {index:.0f} found at {data[index]:.4f}, \t z: {zs[index]:.1f}, \t global p = {p:.5f}, \t reject = {not accept}")
            if accept:
                break
            else:
                mask[index] = False
                mean_copy, std_copy = np.mean(data_copy[mask]), np.std(data_copy[mask], ddof=1)
        
        # Set the mean and std to the new ones
        mean, std = mean_copy, std_copy         
    
    # Plot the new version after applying chauvent's
    ax1.plot(np.arange(len(data)), data, 'r.')
    ax1.errorbar(np.arange(len(data)), data, std, elinewidth = 1, capsize = 2, ls = 'none', c = 'k')
    
    ax1.hlines(mean, -0.5, len(data) + 0.5, label = r"$\mu$")
    ax1.hlines((mean - std, mean + std), -0.5, len(data) + 0.5, ls = '--', label = "$\\mu \\pm \\sigma$")
    
    ax1.legend(loc = "upper left")
    
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
     
    ax.plot(xs, norm.pdf(xs, mean, std), 'k')
    ax.scatter(data, norm.pdf(data, mean, std), c = 'g', zorder = 5, label = "Accepted points")
    ax.scatter(data[np.logical_not(mask)], norm.pdf(data[np.logical_not(mask)], mean, std), c = 'r', label = "Rejected points", zorder = 10)
    
    ax.legend(loc = "upper right")
    
    return fig, [ax, ax1], data_copy[mask]
    
def normalize_pdf(expr, start, end, variable = None):
    """
    Takes sympy expression and return the factor to multiply with to normalize_pdf the function
    """
    from sympy import Symbol, solve, Eq
    C = Symbol("C")

    # If only one variable we just use that
    if variable == None and len(expr.free_symbols) == 1:
        variable = expr.free_symbols

    expr = C * expr

    # Integrate and solve for the factor
    integral = expr.integrate((variable, start, end))
    C = solve(Eq(integral, 1))

    return C
    
def accept_reject(expr, start, end, variable = None):
    """
    Takes a sympy expression, normalize_pdf it and return a function to generate by accept/reject.
    """

    from scipy.optimize import minimize
    from sympy import lambdify

    # Normalize pdf. Not neccessary, but pretty fun
    expr = normalize_pdf(expr, start, end, variable)[0] * expr

    if variable == None and len(expr.free_symbols) == 1:
        x = expr.free_symbols

    lambd = lambdify(x, expr)

    # Find maximum
    max_bound_arg = minimize(lambda a: - lambd(a), x0 = [.5 * (start + end)], bounds = [(start, end)]) 
    
    # Break if no maximum is found
    if not max_bound_arg.success:
        raise RuntimeError("No maximum found")

    
    max_bound = lambd(max_bound_arg.x)

    if max_bound <= 0:
        raise ValueError("Function has negative values")
    
    # generator function. dummy variable is for vectorizing
    def generator(dummy = None):
        while True:
            x = np.random.uniform(start, end)
            y = np.random.uniform(0, max_bound * 1.01) # The multiplication is to make sure we don't undersample
            if y <= lambd(x):
                return x
    
    # Vectorize
    generator_vec = np.vectorize(generator)
    
    # Make it to take a size.
    generator_vec_size = lambda size: generator_vec(np.empty(size))

    return generator_vec_size

from sympy.matrices import Matrix, diag, matrix2numpy
from sympy import Symbol, simplify
class evaluater():
    """
    Klasse til at evaluere en expression eller flere expressions.
    Denne burde indeholde error propagation
    """
    expr = None
    multiple_expr = False
    numerical = False
    variables = []
    covariance = None
    
    def __init__(self, exprs, numerical = False, variables = None, covariance = None, vals = None, errors = None, allow_complex = False):
        # Set expressions. Check if multiple varialbes for future functions
        if type(exprs) not in (list, np.ndarray, tuple):
            self.expr = [exprs]
        else:
            self.expr = exprs

        if len(self.expr) > 1:
            self.multiple_expr = True

        # Check if the results should be numerical
        self.numerical = numerical
        
        # Generate list of variables
        if not variables:
            variables = set([])
            for expr in exprs:
                variables = variables.union(expr.free_symbols)
        self.variables = list(variables)
        
        # make numbers real
        if not allow_complex:
            for var in self.variables:
                var.assumptions0['real'] = True
        
        # Generate covariance matrix
        if covariance:
            self.numerical = True
            self.covariance = Matrix(covariance)
        elif errors:
            self.covariance = np.diag([err ** 2 for err in errors])
        elif variables:
            self.covariance = diag(*[Symbol("\\sigma_{}".format(str(var)), positive = True) ** 2 for var in variables])
        else:
            self.covariance = diag(*[Symbol("\\sigma_{}".format(str(var)), positive = True) ** 2 for var in variables])
        
        if vals:
            self.vals = vals
        
        if vals and errors:
            self.vals = vals
            self.numerical = True
        elif errors and not vals:
            raise InputError("Errors is meaningless without values")
        
    
    def Jacobian(self):
        """ Jacobian matrix. G_ij = df_i/dx_j"""
        print(self.expr)
        return Matrix(len(self.expr), len(self.variables), lambda i, j: self.expr[i].diff(self.variables[j]))
    
    def error_prop_symbolic(self):
        """
        Define Jacobian and calculate the covariance matrix of the transformation. If only one function is given
        the symbolic standard deviation is returned 
        """
        G = self.Jacobian()
        cov = G * self.covariance * G.T
        if not self.multiple_expr:
            return sp.simplify(sp.sqrt(cov[0]))
        return cov.applyfunc(sp.sqrt).applyfunc(simplify)
    
    def error_prop(self):
        """
        The general error_prop function
        """
        if not self.numerical:
            return self.error_prop_symbolic()
        if self.numerical:
            G = self.Jacobian()
            G_np = matrix2numpy(G.subs({i:j for i, j in zip(self.variables, self.vals)}))
            covariance =  (G_np @ self.covariance @ G_np.T)
            if not self.multiple_expr:
                return np.sqrt(float(covariance[0][0]))

    
if __name__ == "__main__":
    from sympy.abc import x
    expr = - 0.25 * x + 15 + 15 * sp.exp(- (x - 55) ** 2 / (2.5) ** 2 / 2)
    xmin, xmax = 0, 100

    gen = accept_reject(expr, xmin, xmax)
    # print("generated function")

    data = gen(10000)
    # print(data)
    fig, ax = plt.subplots()
    ax.hist(data, bins = 50, histtype = 'step')

    fig.show()

    # print(accept_reject(expr, xmin, xmax))
