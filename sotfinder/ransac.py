import numpy as np
from sklearn.metrics import r2_score
from sklearn.linear_model import RANSACRegressor

class PolynomialRegression(object):
    """
    Polynomial Regression model for fitting a quadratic curve to data.
    """

    def __init__(self, degree=2, coeffs=None):
        """
        Initializes the PolynomialRegression model with default parameters.
        """
        self.degree = degree
        self.coeffs = coeffs

    def fit(self, X, y):
        """
        Fits a quadratic polynomial to the provided data.

        Inputs:
            X -> [array-like] Independent variable.
            y -> [array-like] Dependent variable.
        """

        # Fit the model coefficients
        self.coeffs = np.polyfit(X.ravel(), y, self.degree)

    def get_params(self, deep=False):
        """
        Returns the coefficients of the polynomial regression model.
        """
        return {'coeffs': self.coeffs}   

    def set_params(self, coeffs=None, random_state=None):
        """
        Sets the coefficients of the polynomial regression model.
        """
        self.coeffs = coeffs     

    def predict(self, X):
        """
        Make prediction using the fitted polynomial model.

        Inputs:
            X -> [array-like] Independent variable for which to make predictions.
        Outputs:
            Predicted values (array-like).
        """
        poly_eqn = np.poly1d(self.coeffs)
        y_hat = poly_eqn(X.ravel())
        return y_hat

    def score(self, X, y):  
        """ 
        Calculates the R² score, a measure of how well the predicted values match the actual values.

        Inputs:
            X -> [array-like] Independent variable.
            y -> [array-like] Dependent variable.
        Outputs:
            score -> [float] The R² score.
        """    
        return r2_score(y,self.predict(X))  
        
def sequential_ransac(x,y,min_inliers_allowed=10,min_samples=3,residual_threshold=5.0):
    """
    Sequentially applies the RANSAC algorithm to identify multiple motion trajectories.

    Inputs:
        X -> [array-like] Independent variables (e.g., x coordinates).
        y -> [array-like] Dependent variables (e.g., y coordinates).
        min_inliers_allowed -> [int,optional,default=10] Minimum number of remaining points required for another RANSAC iteration.
        min_samples -> [int,optional,default=3] Minimum number of data points required to fit the model.
        residual_multiplier -> [float,optional,default=5.0] Multiplier for MAD(Median Absolute Deviation) sigma to set the residual threshold.
    Outputs:
        lines -> [dict] A dictionary representing detected trajectories.
    """
    # Initialize the RANSAC regressor with the given model
    ransac_model = PolynomialRegression()
    ransac = RANSACRegressor(ransac_model,min_samples=min_samples,residual_threshold=residual_threshold,max_trials=1000)

    X = x[:,None]

    lines = {} # Initialize a list to store the detected trajectories
    i = 0 # Sequence of trajectories
    # Continue until there are not enough points left for next iteration
    while len(X) > min_inliers_allowed:
        ransac.fit(X, y)
        # Determine inliers and outliers
        inlier_mask = ransac.inlier_mask_
        outlier_mask = ~inlier_mask 

        # Break if no inliers to form a new trajectory
        if not inlier_mask.any(): 
            break 
        
        # Update data with remaining outliers for the next iteration
        X_inlier,X_outlier = X[inlier_mask],X[outlier_mask]
        y_inlier,y_outlier = y[inlier_mask],y[outlier_mask]
        lines[str(i)] = np.stack([X_inlier.ravel(),y_inlier]).T
        X,y = X_outlier,y_outlier
        i += 1
    return lines # Return the list of detected trajectories