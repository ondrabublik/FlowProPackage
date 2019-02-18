/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.user.domaintransformationobject;

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
}
