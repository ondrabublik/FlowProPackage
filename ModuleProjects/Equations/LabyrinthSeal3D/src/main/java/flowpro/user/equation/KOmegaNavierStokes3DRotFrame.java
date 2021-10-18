/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.user.equation;

import flowpro.api.ElementData;
import flowpro.api.FlowProProperties;
import flowpro.api.Mat;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class KOmegaNavierStokes3DRotFrame extends NavierStokes3DRotFrame {
    
    protected static final double P_TOL = 1e-1;

    double kIn;
    double omIn;
    double omWall;
    
    double kref;
    double omref;
    double muref;

    double LIM_TOL;

    // Sutherland constant
    double Suth;
    
    // parameters of turbulence model
    double a1;
    double sk1;
    double som1;
    double alpha1;
    double beta1;
    double betast;   // beta star 
    double sk2;
    double som2;
    double alpha2;
    double beta2;
    double Prt;

    public void init(FlowProProperties props) throws IOException {
        int dimension = props.getInt("dimension");
        super.init(props);
        this.nEqs = 6;

        // inlet
        kIn = props.getDouble("kIn");
        omIn = props.getDouble("omIn");
        if (isInletSupersonic) {
            WIn[dimension + 2] = kIn;
            WIn[dimension + 3] = omIn;
        }
        
        // Sutherland
        Suth = props.getDouble("SuthConst");
        
        // parameters of turbulence model
        a1 = props.getDouble("a1");
        sk1 = props.getDouble("sk1");
        som1 = props.getDouble("som1");
        alpha1 = props.getDouble("alpha1");
        beta1 = props.getDouble("beta1");
        betast = props.getDouble("betast");
        sk2 = props.getDouble("sk2");
        som2 = props.getDouble("som2");
        alpha2 = props.getDouble("alpha2");
        beta2 = props.getDouble("beta2");
        Prt = props.getDouble("Prt");
    }

    @Override
    public double[] constInitCondition() {
        if (isInletSupersonic) {
            return WIn;
        } else {
            double machIn = Math.sqrt(2 / (kapa - 1) * (Math.pow((1 / pOut), (kapa - 1) / kapa) - 1));
            double rhoIn = Math.pow(1 + ((kapa - 1) / 2) * machIn * machIn, 1 / (1 - kapa));
            double VIn = machIn * Math.sqrt((kapa * pOut) / rhoIn);
            double uIn = VIn * Math.cos(attackAngleAlfa) * Math.cos(attackAngleBeta);
            double vIn = VIn * Math.sin(attackAngleAlfa) * Math.cos(attackAngleBeta);
            double wIn = VIn * Math.sin(attackAngleAlfa) * Math.sin(attackAngleBeta);
            double EIn = pOut / (kapa - 1) + 0.5 * rhoIn * VIn * VIn;

            return new double[]{rhoIn, rhoIn * uIn, rhoIn * vIn, rhoIn * wIn, EIn, rhoIn * kIn, rhoIn * omIn};
        }
    }

    //  nevazky tok stenou _____________________________________________________
    @Override
    public double[] numericalConvectiveFlux(double WL[], double WR[], double[] n, int TT, ElementData elem) {
        WL[0] = limiteRho(WL[0]);
        WR[0] = limiteRho(WR[0]);

        double[] fn = new double[nEqs];
        double p;

        switch (TT) {
            case (NavierStokes3DRotFrame.BoundaryType.WALL):
            case (NavierStokes3DRotFrame.BoundaryType.INVISCID_WALL):
            case (NavierStokes3DRotFrame.BoundaryType.STATOR):
            case (NavierStokes3DRotFrame.BoundaryType.ROTOR): // stena
                p = pressure(WL, elem.currentX);
                fn[0] = 0;
                fn[1] = p * n[0];
                fn[2] = p * n[1];
                fn[3] = p * n[2];
                fn[4] = 0;
                fn[5] = 0;
                fn[6] = 0;
                break;

            case (NavierStokes3DRotFrame.BoundaryType.INLET): // vstup
                fn = convectiveFlux(WR, n, elem);
                break;

            case (NavierStokes3DRotFrame.BoundaryType.OUTLET): // vystup
                fn = convectiveFlux(WR, n, elem);
                break;

            default: // vnitrni stena
                double pL = pressure(WL, elem.currentX);
                double pR = pressure(WR, elem.currentX);
                double aL = Math.sqrt(kapa * pL / WL[0]);
                double aR = Math.sqrt(kapa * pR / WR[0]);
                double VnL = WL[1] / WL[0] * n[0] + WL[2] / WL[0] * n[1] + WL[3] / WL[0] * n[2];
                double VnR = WR[1] / WR[0] * n[0] + WR[2] / WR[0] * n[1] + WR[3] / WR[0] * n[2];
                double[] fcL = convectiveFlux(WL, n, elem);
                double[] fcR = convectiveFlux(WR, n, elem);
                double S = Math.max(Math.abs(VnL) + aL, Math.abs(VnR) + aR);
                for (int j = 0; j < nEqs; j++) {
                    fn[j] = (fcL[j] + fcR[j] - S * (WR[j] - WL[j])) / 2;
                }
                break;
        }
        return fn;
    }

    @Override
    public double[] convectiveFlux(double[] W, double[] n, ElementData elem) {
        W[0] = limiteRho(W[0]);

        double[] f = new double[nEqs];
        double V = W[1] / W[0] * n[0] + W[2] / W[0] * n[1] + W[3] / W[0] * n[2];
        double p = pressure(W, elem.currentX);
        f[0] = W[0] * V;
        f[1] = W[1] * V + p * n[0];
        f[2] = W[2] * V + p * n[1];
        f[3] = W[3] * V + p * n[2];
        f[4] = (W[4] + p) * V;
        f[5] = W[5] * V;
        f[6] = W[6] * V;

        return f;
    }

    @Override
    public double[] diffusiveFlux(double[] W, double[] dW, double[] n, ElementData elem) {
        W[0] = limiteRho(W[0]);
        double rho = W[0];

        double lam = -2. / 3; // Stokesuv vztah

        double[] velocity = new double[dim];
        double velocity2 = .0;
        for (int d = 0; d < dim; ++d) {
            velocity[d] = W[d + 1] / rho;
            velocity2 += velocity[d] * velocity[d];
        }

        double[] velocityJac = new double[dim * dim];
        for (int d = 0; d < dim; ++d) {
            for (int f = 0; f < dim; ++f) {
                velocityJac[dim * d + f] = (dW[f * nEqs + d + 1] - dW[f * nEqs] * velocity[d]) / rho;
            }
        }

        double p = pressure(W, elem.currentX);
        double[] pOverRhoDer = new double[dim];
        for (int d = 0; d < dim; ++d) {
            double temp = .0;
            for (int f = 0; f < dim; ++f) {
                temp += velocity[f] * velocityJac[dim * f + d];
            }
            double pDer = (kapa - 1) * (dW[d * nEqs + dim + 1] - dW[d * nEqs] * velocity2 / 2 - rho * temp);
            pOverRhoDer[d] = (rho * pDer - p * dW[d * nEqs]) / (rho * rho);
        }

        // stress tensor calculation
        double[] stress = new double[dim * dim];
        double trace = .0;
        for (int d = 0; d < dim; ++d) {
            trace += velocityJac[dim * d + d];
            for (int f = 0; f < dim; ++f) {
                stress[dim * d + f] = velocityJac[dim * d + f] + velocityJac[dim * f + d];
            }
        }
        for (int d = 0; d < dim; ++d) {
            stress[dim * d + d] += lam * trace;
        }

        double k = max(W[dim + 2]/rho,1e-10);
        double om = W[dim + 3] / rho;
        
        // sutherland relation
        double etaTemp = sutherland(rho, p);
        
        double vortic = Math.abs(velocityJac[1]-velocityJac[2]);
        double walldist = Math.abs(elem.currentWallDistance);  
        if (walldist < 1e-5){                // zero at wall makes problems - adjust value according to mesh
            walldist = 1e-5;
        }
        double arg2 = max(2*Math.sqrt(k)/(0.09*Math.exp(om)*walldist)/Math.sqrt(Re),500*etaTemp/rho/(walldist*walldist*Math.exp(om))/Re);
        double F2 = Math.tanh(arg2*arg2);
        
        double mut = max(0, rho*a1*k / max(a1*Math.exp(om),vortic*F2));

        double[] turbulentStress = new double[dim * dim];
        for (int i = 0; i < stress.length; i++) {
            turbulentStress[i] = stress[i] * mut;
        }
        double kMax = 2.0 / 3 * max(0, rho*k);
        for (int d = 0; d < dim; ++d) {
            turbulentStress[dim * d + d] -= kMax;
        }

        double[] kDer = new double[dim];
        double[] omDer = new double[dim];
        for (int d = 0; d < dim; ++d) {
            kDer[d] = (dW[d * nEqs + dim + 2] - dW[d * nEqs] * k) / rho;
            omDer[d] = (dW[d * nEqs + dim + 3] - dW[d * nEqs] * om) / rho;
        }

        double constant = kapa / (kapa - 1) * (etaTemp / Pr + mut / Prt);
        
        double omkDer = 0;
        for (int d = 0; d < dim; ++d) {
            omkDer += omDer[d] * kDer[d];
        }
        
        double CD = max(2*som2*rho*omkDer,1e-10);
        double arg11 = max(Math.sqrt(k)/(0.09*Math.exp(om)*walldist)/Math.sqrt(Re),500*etaTemp/rho/(Math.exp(om)*walldist*walldist)/Re);
        double arg12 = 4*rho*som2*k/(CD*walldist*walldist);
        double arg1 = min(arg11,arg12);   
        double F1 = Math.tanh(Math.pow(arg1,4));
        
        double sk0 = F1*sk1 + (1-F1)*sk2;
        double som0 = F1*som1 + (1-F1)*som2;
        
        double[] flux = new double[nEqs];
        for (int d = 0; d < dim; ++d) {
            double tmp = .0;
            for (int f = 0; f < dim; ++f) {
                flux[f + 1] += (etaTemp * stress[dim * d + f] + turbulentStress[dim * d + f]) * n[d] / Re;
                tmp += velocity[f] * (etaTemp * stress[dim * d + f] + turbulentStress[dim * d + f]);
            }
            flux[dim + 1] += (tmp + constant * pOverRhoDer[d]) * n[d] / Re;
            flux[dim + 2] += (etaTemp + sk0 * mut) / Re * kDer[d] * n[d];
            flux[dim + 3] += (etaTemp + som0 * mut) / Re * omDer[d] * n[d];
        }
        
        return flux;
    }

    @Override
    public double[] sourceTerm(double[] W, double[] dW, ElementData elem) { // zdrojovy clen
        W[0] = limiteRho(W[0]);

        double[] source = new double[nEqs];
        double[] Fcentrifugal = Mat.times(Mat.cross(omega, Mat.cross(omega, elem.currentX)), -W[0]);
        double[] Wref = new double[]{W[1] / W[0], W[2] / W[0], W[3] / W[0]};
        double[] Fcoriolis = Mat.times(Mat.cross(omega, Wref), -2.0 * W[0]);
        source[1] = Fcentrifugal[0] + Fcoriolis[0];
        source[2] = Fcentrifugal[1] + Fcoriolis[1];
        source[3] = Fcentrifugal[2] + Fcoriolis[2];

        W[0] = limiteRho(W[0]);
        double rho = W[0];

        double lam = -2. / 3; // Stokesuv vztah

        double[] velocity = new double[dim];
        for (int d = 0; d < dim; ++d) {
            velocity[d] = W[d + 1] / rho;
        }
        
        double[] velocityJac = new double[dim * dim];
        for (int d = 0; d < dim; ++d) {
            for (int f = 0; f < dim; ++f) {
                velocityJac[dim * d + f] = (dW[f * nEqs + d + 1] - dW[f * nEqs] * velocity[d]) / rho;
            }
        }
        
        double p = pressure(W);

        double k = max(W[dim + 2]/rho,1e-10);
        double om = W[dim + 3] / rho;
        double expOm = Math.exp(om);
        
        double etaTemp = sutherland(rho, p);
        
        double vortic = Math.abs(velocityJac[1]-velocityJac[2]);
        double walldist = Math.abs(elem.currentWallDistance);  
        if (walldist < 1e-5){                // zero at wall makes problems - adjust value according to mesh
            walldist = 1e-5;
        }
        double arg2 = max(2*Math.sqrt(k)/(0.09*Math.exp(om)*walldist)/Math.sqrt(Re),500*etaTemp/rho/(walldist*walldist*Math.exp(om))/Re);
        double F2 = Math.tanh(arg2*arg2);
        
        double mut = max(0, rho*a1*k / max(a1*Math.exp(om),vortic*F2));      

        // turbulentStress tensor calculation
        double[] turbulentStress = new double[dim * dim];
        double trace = .0;
        for (int d = 0; d < dim; ++d) {
            trace += velocityJac[dim * d + d];
            for (int f = 0; f < dim; ++f) {
                turbulentStress[dim * d + f] = mut * (velocityJac[dim * d + f] + velocityJac[dim * f + d]);
            }
        }
        double kMax = 2.0 / 3 * max(0, rho*k);
        for (int d = 0; d < dim; ++d) {
            turbulentStress[dim * d + d] += mut * lam * trace - kMax;
        }

        double omkDer = 0;
        double omDerSqr = 0;
        for (int d = 0; d < dim; ++d) {
            double kDer = (dW[d * nEqs + dim + 2] - dW[d * nEqs] * k) / rho;
            double omDer = (dW[d * nEqs + dim + 3] - dW[d * nEqs] * om) / rho;
            omkDer += omDer * kDer;
            omDerSqr += omDer * omDer;
        }

        double Tv = 0;
        for (int d = 0; d < dim; d++) {
            for (int f = 0; f < dim; f++) {
                Tv += turbulentStress[dim * d + f] * velocityJac[dim * d + f];
            }
        }
              
        double CD = max(2*som2*rho*omkDer,1e-10);
        double arg11 = max(Math.sqrt(k)/(0.09*Math.exp(om)*walldist)/Math.sqrt(Re),500*etaTemp/rho/(Math.exp(om)*walldist*walldist)/Re);
        double arg12 = 4*rho*som2*k/(CD*walldist*walldist);
        double arg1 = min(arg11,arg12);   
        double F1 = Math.tanh(Math.pow(arg1,4));
        
        //double sk0 = F1*sk1 + (1-F1)*sk2;
        double som0 = F1*som1 + (1-F1)*som2;
        double alpha0 = F1*alpha1 + (1-F1)*alpha2;
        double beta0 = F1*beta1 + (1-F1)*beta2;

        if (Tv > 10 * betast * rho * k * expOm) {   
            Tv = 10 * betast * rho * k * expOm;
        }

        double Tvom = alpha0 * rho / mut /expOm * Tv;   
        if (Tvom > 10 * beta0 * rho * expOm) {     // moje zmena, opatrne (bylo 100*...)
            Tvom = 10 * beta0 * rho * expOm;
        }

        source[dim + 2] = max(Tv - betast * rho * k * expOm, 0);  
        source[dim + 3] = max(Tvom - beta0 * rho * expOm + 1 / Re * som2*rho/expOm *omkDer *2*(1-F1) + 1 / Re * (etaTemp + som0 * mut) * omDerSqr, 0);
        return source;
    }

    @Override
    public double[] boundaryValue(double[] WL, double[] n, int TT, ElementData elem) {
        WL[0] = limiteRho(WL[0]);

        double[] WR = new double[nEqs];
        double p = pressure(WL, elem.currentX);
        switch (TT) {
            case (NavierStokes3DRotFrame.BoundaryType.WALL): // stena
                if (isDiffusive) {
                    WR[0] = WL[0];
                    WR[1] = 0;
                    WR[2] = 0;
                    WR[3] = 0;
                    WR[4] = p / (kapa - 1);
                    WR[5] = 0;
                    WR[6] = WL[dim + 3];
                } else {
                    System.arraycopy(WL, 0, WR, 0, nEqs);
                }
                break;

            case (NavierStokes3DRotFrame.BoundaryType.STATOR): // stator
                if (isDiffusive) {
                    double[] uRot = Mat.cross(omegaStator, elem.currentX);
                    WR[0] = WL[0];
                    WR[1] = WR[0] * uRot[0];
                    WR[2] = WR[0] * uRot[1];
                    WR[3] = WR[0] * uRot[2];
                    WR[4] = p / (kapa - 1) + (WR[1] * WR[1] + WR[2] * WR[2] + WR[3] * WR[3]) / (2 * WR[0]) - WR[0] * Math.pow(Mat.L2Norm(Mat.cross(omega, elem.currentX)), 2) / 2;
                    WR[5] = 0;
                    WR[6] = WL[dim + 3];
                } else {
                    System.arraycopy(WL, 0, WR, 0, nEqs);
                }
                break;

            case (NavierStokes3DRotFrame.BoundaryType.ROTOR): // rotor
                if (isDiffusive) {
                    double[] uRot = Mat.cross(omegaRotor, Mat.minusVec(elem.currentX, rotorExcentricity));
                    WR[0] = WL[0];
                    WR[1] = WR[0] * uRot[0];
                    WR[2] = WR[0] * uRot[1];
                    WR[3] = WR[0] * uRot[2];
                    WR[4] = p / (kapa - 1) + (WR[1] * WR[1] + WR[2] * WR[2] + WR[3] * WR[3]) / (2 * WR[0]) - WR[0] * Math.pow(Mat.L2Norm(Mat.cross(omega, elem.currentX)), 2) / 2;
                    WR[5] = 0;
                    WR[6] = WL[dim + 3];
                } else {
                    System.arraycopy(WL, 0, WR, 0, nEqs);
                }
                break;

            case (NavierStokes3DRotFrame.BoundaryType.INLET):
                if (WIn[0] == -1) { // subsonicky vstup
                    double Minl;
                    if (p > pIn0) {
                        Minl = 0;
                    } else {
                        Minl = Math.sqrt((2 / (kapa - 1)) * (-1 + Math.pow(pIn0 / p, (kapa - 1) / kapa)));
                    }
                    double[] uRot = new double[3];
                    if (swirlProfileRadius != null) {
                        double[] uRotR = Mat.cross(omegaRotor, elem.currentX);
                        double[] uRotS = Mat.cross(omegaStator, elem.currentX);
                        double[] uRotIn = Mat.cross(omegaInlet, elem.currentX);
                        double rR = swirlProfileRadius[0]; // rotor radius
                        double rS = swirlProfileRadius[1]; // stator radius
                        double rIn = (rR + rS) / 2;
                        double x = Mat.L2Norm(elem.currentX);
                        for (int i = 0; i < 3; i++) {
                            uRot[i] = uRotR[i]*(x-rIn)*(x-rS)/(rR-rIn)/(rR-rS) + uRotIn[i]*(x-rR)*(x-rS)/(rIn-rR)/(rIn-rS) + uRotS[i]*(x-rR)*(x-rIn)/(rS-rR)/(rS-rIn);
                        }
                    } else {
                        uRot = Mat.cross(omegaInlet, elem.currentX);
                    }
                    double Rinl = rhoIn0 * Math.pow((1 + ((kapa - 1) / 2) * Minl * Minl), 1 / (1 - kapa));
                    double Vinl = Minl * Math.sqrt((kapa * p) / Rinl);
                    double uinl = Vinl * Math.cos(attackAngleAlfa) * Math.cos(attackAngleBeta) + uRot[0];
                    double vinl = Vinl * Math.sin(attackAngleAlfa) * Math.cos(attackAngleBeta) + uRot[1];
                    double winl = Vinl * Math.sin(attackAngleAlfa) * Math.sin(attackAngleBeta) + uRot[2];
                    double Einl = p / (kapa - 1) + 0.5 * Rinl * Vinl * Vinl - 0*Rinl * Math.pow(Mat.L2Norm(Mat.cross(omega, elem.currentX)), 2) / 2;
                    WR[0] = Rinl;
                    WR[1] = Rinl * uinl;
                    WR[2] = Rinl * vinl;
                    WR[3] = Rinl * winl;
                    WR[4] = Einl;
                    WR[5] = Rinl * kIn;
                    WR[6] = Rinl * omIn;
                } else { // supersonicky vstup
                    WR[0] = WIn[0];
                    WR[1] = WIn[1];
                    WR[2] = WIn[2];
                    WR[3] = WIn[3];
                    WR[4] = WIn[4];
                    WR[5] = WIn[0] * WIn[5];
                    WR[6] = WIn[0] * WIn[6];
                }
                break;

            case (NavierStokes3DRotFrame.BoundaryType.OUTLET):
                double ro = WL[0];
                double uo = WL[1] / WL[0];
                double vo = WL[2] / WL[0];
                double wo = WL[3] / WL[0];
                double kOut = WL[5] / WL[0];
                double omOut = WL[5] / WL[0];
                double ao = Math.sqrt(kapa * p / ro);
                double Mo = Math.sqrt(uo * uo + vo * vo + wo * wo) / ao;
                double Eo,
                 pout;
                if (Mo < 1) { // subsonicky vystup
                    pout = pOut;
                } else { // supersonicky vystup
                    pout = p;
                }
                Eo = pout / (kapa - 1) + ro * (uo * uo + vo * vo + wo * wo) / 2 - ro * Math.pow(Mat.L2Norm(Mat.cross(omega, elem.currentX)), 2) / 2;
                WR[0] = ro;
                WR[1] = ro * uo;
                WR[2] = ro * vo;
                WR[3] = ro * wo;
                WR[4] = Eo;
                WR[5] = ro * kOut;
                WR[5] = ro * omOut;
                break;

            case (NavierStokes3DRotFrame.BoundaryType.INVISCID_WALL): // nevazka stena
                System.arraycopy(WL, 0, WR, 0, nEqs);
                break;
        }
        return WR;
    }

    public double sutherland(double rho, double p) {
        double T = p / rho;
        double TRef = 1/ (cv*(kapa - 1)) * pRef / rhoRef;
        return Math.pow(T,1.5)*(1+Suth/TRef)/(T + Suth/TRef);
    }
    
    public double max(double a, double b) {
        if (a > b) {
            return a;
        } else {
            return b;
        }
    }
    
    public double min(double a, double b) {
        if (a < b) {
            return a;
        } else {
            return b;
        }
    }

    public double matrixMagnitude(double[] A) {
        double mag = 0;
        for (int i = 0; i < A.length; i++) {
            mag += A[i] * A[i];
        }
        return Math.sqrt(mag);
    }

    public double rotationMagnitude(double[] U) {
        double rotMag = 0;
        if (U.length == 4) {
            rotMag = Math.abs(U[1] - U[2]);
        } else if (U.length == 9) {
            double w1 = (U[0*dim + 1] - U[1*dim + 0]);
            double w2 = (U[0*dim + 2] - U[2*dim + 0]);
            double w3 = (U[2*dim + 1] - U[1*dim + 2]);
            rotMag = Math.sqrt(w1*w1 + w2*w2 + w3*w3);
        }
        return rotMag;
    }
    
    public double limitDestruction(double d) {
        if (d < -100) {
            d = -100;
        }
        return d;
    }
    
    @Override
    public boolean isSourcePresent() {
        return true;
    }

    @Override
    public boolean isIPFace(int TT) {
        if (TT == NavierStokes3DRotFrame.BoundaryType.WALL || TT == NavierStokes3DRotFrame.BoundaryType.STATOR || TT == NavierStokes3DRotFrame.BoundaryType.ROTOR) {
            return true;
        } else {
            return false;
        }
    }
    
    @Override
    public double[] getResults(double[] W, double[] dW, double[] X, String name) {
        double[] uRot;
        switch (name) {
            case "mach":
                double absVelocity = .0;
                for (int d = 0; d < dim; ++d) {
                    absVelocity += W[d + 1] * W[d + 1];
                }
                absVelocity = Math.sqrt(absVelocity) / W[0];

                double a = Math.sqrt(kapa * pressure(W) / W[0]);
                return new double[]{absVelocity / a};

            case "density":
                return new double[]{rhoRef * W[0]};

            case "xVelocity":
                uRot = Mat.cross(omega, X);
                return new double[]{velocityRef * (W[1] / W[0] + uRot[0])};

            case "yVelocity":
                uRot = Mat.cross(omega, X);
                return new double[]{velocityRef * (W[2] / W[0] + uRot[1])};

            case "zVelocity":
                uRot = Mat.cross(omega, X);
                return new double[]{velocityRef * (W[3] / W[0] + uRot[2])};

            case "velocity":
                uRot = Mat.cross(omega, X);
                double[] velocity = new double[dim];
                for (int i = 0; i < dim; i++) {
                    velocity[i] = velocityRef * (W[i + 1] / W[0] + uRot[i]);
                }
                return velocity;

            case "temperature":
                throw new UnsupportedOperationException("undefined value" + name);

            case "energy":
                return new double[]{pRef * W[dim + 1]};

            case "pressure":
                return new double[]{pRef * pressure(W)};
                
            case "mut":
                return new double[]{W[dim + 2] / Math.exp(W[dim + 3] / W[0])};
                
            default:
                throw new UnsupportedOperationException("undefined value " + name);
        }
    }
}
