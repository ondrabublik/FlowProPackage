/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.user.domaintransformationobject;

import flowpro.api.Complex;
import static flowpro.api.Complex.*;
import flowpro.api.DomainTransformationObject;
import flowpro.api.FlowProProperties;

/**
 *
 * @author obublik
 */
public class CylindricalTransform implements DomainTransformationObject{
    
    public void init(FlowProProperties props){
        
    }
    
    public double[] transform(double[] X) {
        return new double[]{X[0], X[2] * Math.cos(X[1]), -X[2] * Math.sin(X[1])};
    }
    
    public Complex[] transformComplex(Complex[] X){
        Complex[] Y = new Complex[X.length];
        Y[0] = copy(X[0]);
        Y[1] = multiply(X[2],cos(X[1]));
        Y[2] = multiply(multiply(X[2],sin(X[1])),new Complex(-1,0));
        
        return Y;
    }
}
