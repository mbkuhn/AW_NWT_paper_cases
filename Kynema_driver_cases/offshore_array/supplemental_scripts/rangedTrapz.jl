using Interpolations
#using BSplineKit
using Trapz
function rangedTrapz(x,y,lim)
        #Function computing the integral of a data sequence over the specified range.
        #   -  Matteo Polimeno, NREL, 15 April 2025
        #   - Adapted from MATLAB's implementation by Lu Wang, NREL, 10/17/2019
        
        # Inputs: x - Independent variable of the data sequence
        #         y - Dependent variable of the data sequence
        #         lim - 2-by-1 array. lim(1) and lim(2) give the lower and upper
        #               bounds of the x range over which the integral is to be computed.
        
        # Outputs: int - Integral of y over the specified x range.
        
        # Check if limit is out of range
        if lim[1]<x[1] || lim[2]>x[end]
            error("rangedMax: Specified limits are out of range") 
        end
        
        # Extract relevant portion of the data.
        # Determine endpoint values using linear interpolation.
        startIdx_all = findall(x->x<=lim[1],x)
        startIdx = startIdx_all[end]
        endIdx_all = findall(x->x>=lim[2],x)
        endIdx = endIdx_all[1]
        x0 = vec(x[startIdx:endIdx])
        y0 = vec(y[startIdx:endIdx])
        x0[1] = lim[1] 
        x0[end] = lim[2]
        y_f = linear_interpolation(x0,y0)
        #y_f = interpolate(x0,y0,BSplineOrder(4))
        yStart = y_f(lim[1])
        yEnd = y_f(lim[2])
        y0[1] = yStart 
        y0[end] = yEnd
        
        #Integration with trapezoidal rule
        int = trapz(x0,y0)
        
        return int
    
    end