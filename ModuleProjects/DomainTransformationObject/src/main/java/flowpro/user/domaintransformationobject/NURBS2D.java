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
public class NURBS2D implements DomainTransformationObject {
    
    public void init(FlowProProperties props){
        
    }
    
    public double[] transform(double[] X){
        double tx = X[0];
        double ty = X[1];
        
        double x = tx;
        double bump = 0.1 * Math.exp(-(tx - 1.5) * (tx - 1.5) / 0.02);
        double y = ty + bump*(1-ty);
        
//        double x = 3*tx;
//        double bump = 0.1 * Math.exp(-(tx - 0.5) * (tx - 0.5) / 0.02);
//        double y = (1-bump)*ty + bump*(1-ty);
        return new double[]{x,y};
        //return X;
    }
}
