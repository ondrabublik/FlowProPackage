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
public class SpalartAllmarasNavierStokes3DRotFrame extends NavierStokes3DRotFrame {
    
    double vtIn; // turbulence intensity at the inlet

    // model turbulence
    int model;
    double sigma = 2.0 / 3;
    double cb1 = 0.1355;
    double cb2 = 0.622;
    double ka = 0.41;
    double cw1 = cb1 / (ka * ka) + (1 + cb2) / sigma;
    double cw2 = 0.3;
    double cw3 = 2;
    double cv1 = 7.1;
    double ct1 = 1;
    double ct2 = 2;
    double ct3 = 1.2;
    double ct4 = 0.5;
    double Prt = 0.9;
    double C_prod = 2;

    public void init(FlowProProperties props) throws IOException {

        // inlet
        vtIn = props.getDouble("vtIn");

        // parameters of turbulence model
        sigma = props.getDouble("sigma");
        cb1 = props.getDouble("cb1");
        cb2 = props.getDouble("cb2");
        ka = props.getDouble("ka");
        cw1 = cb1 / (ka * ka) + (1 + cb2) / sigma;
        cw2 = props.getDouble("cw2");
        cw3 = props.getDouble("cw3");
        cv1 = props.getDouble("cv1");
        ct1 = props.getDouble("ct1");
        ct2 = props.getDouble("ct2");
        ct3 = props.getDouble("ct3");
        ct4 = props.getDouble("ct4");
        Prt = props.getDouble("Prt");
        C_prod = props.getDouble("C_prod");

        super.init(props);
        this.nEqs = 6;
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

            return new double[]{rhoIn, rhoIn * uIn, rhoIn * vIn, rhoIn * wIn, EIn, rhoIn * vtIn};
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
            case (BoundaryType.WALL):
            case (BoundaryType.INVISCID_WALL):
            case (BoundaryType.STATOR):
            case (BoundaryType.ROTOR): // stena
                p = pressure(WL, elem.currentX);
                fn[0] = 0;
                fn[1] = p * n[0];
                fn[2] = p * n[1];
                fn[3] = p * n[2];
                fn[4] = 0;
                fn[5] = 0;
                break;

            case (BoundaryType.INLET): // vstup
                fn = convectiveFlux(WR, n, elem);
                break;

            case (BoundaryType.OUTLET): // vystup
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

        double vt = W[5] / rho;
        if(vt < 0){
            vt = 0;
        }
        double[] dvt = new double[dim]; 
        for (int d = 0; d < dim; ++d) {
            dvt[d] = (dW[5] - dW[d * nEqs] * vt)/rho;
        }
        
        double xi = rho * vt; // vt/v
        if(xi > 1e5){
            xi = 1e5;
        }
        double fv1 = xi * xi * xi / (xi * xi * xi + cv1 * cv1 * cv1);
        double mut = rho * vt * fv1;
        
        // add turbulence stress
        for (int d = 0; d < dim; ++d) {
            for (int f = 0; f < dim; ++f) {
                stress[dim * d + f] *= (1 + mut);
            }
        }
        
        double constant = kapa / (kapa - 1)*(1/ Pr + mut / Prt);
        double[] flux = new double[nEqs];
        flux[0] = 0;
        for (int d = 0; d < dim; ++d) {
            double tmp = .0;
            for (int f = 0; f < dim; ++f) {
                flux[f + 1] += stress[dim * d + f] * n[d] / Re;
                tmp += velocity[f] * stress[dim * d + f];
            }
            flux[dim + 1] += (tmp + constant * pOverRhoDer[d]) * n[d] / Re;
            flux[dim + 2] += 1 / (Re * sigma) * (1 + rho * vt) * dvt[d] * n[d];
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

        // sa model
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
        double Smag = matrixMagnitude(stress) / 2;

        double vt = max(0, W[dim + 2] / rho);
        double vtDerMag = 0;
        for (int d = 0; d < dim; ++d) {
            double vtDer = 1 / rho * (dW[d * nEqs + dim + 2] - dW[d * nEqs] * vt);
            vtDerMag += vtDer * vtDer;
        }

        // turbulence limit
        if (vt < 0) {
            vt = 0;
        }

        double D = elem.currentWallDistance;

        double xi = rho*vt; // vt/v
        if(xi > 1e5){
            xi = 1e5;
        }
        double ft2 = 0; //ct3*Math.exp(-ct4*xi*xi);
        double fv1 = xi * xi * xi / (xi * xi * xi + cv1 * cv1 * cv1);
        double fv2 = 1 - xi / (1 + xi * fv1);
        double Om = rotationMagnitude(velocityJac);
        // S = Om + C_prod * Math.min(0, Smag - Om) + 1 / Re * vt / (ka * ka * D * D) * fv2;
        double S = Om + 1 / Re * vt / (ka * ka * D * D) * fv2;
        double rt = vt / (Re * S * ka * ka * D * D);
        if (rt > 10) {
            rt = 10;
        }
        double g = rt + cw2 * (Math.pow(rt, 6.0) - rt);
        double fw = g * Math.pow((1 + Math.pow(cw3, 6.0)) / (Math.pow(g, 6.0) + Math.pow(cw3, 6.0)), 1.0 / 6);
        
        source[dim + 2] = limitDestruction(1 / Re * rho * cb2 * vtDerMag + cb1 * (1 - ft2) * rho * S * vt - 1 / Re * (cw1 * fw - cb1 / (ka * ka) * ft2) / rho * (rho * vt / D) * (rho * vt / D));
        return source;
    }

    @Override
    public double[] boundaryValue(double[] WL, double[] n, int TT, ElementData elem) {
        WL[0] = limiteRho(WL[0]);

        double[] WR = new double[nEqs];
        double p = pressure(WL, elem.currentX);
        switch (TT) {
            case (BoundaryType.WALL): // stena
                if (isDiffusive) {
                    WR[0] = WL[0];
                    WR[1] = 0;
                    WR[2] = 0;
                    WR[3] = 0;
                    WR[4] = p / (kapa - 1);
                    WR[5] = 0;
                } else {
                    System.arraycopy(WL, 0, WR, 0, nEqs);
                }
                break;

            case (BoundaryType.STATOR): // stator
                if (isDiffusive) {
                    double[] uRot = Mat.cross(omegaStator, elem.currentX);
                    WR[0] = WL[0];
                    WR[1] = WR[0] * uRot[0];
                    WR[2] = WR[0] * uRot[1];
                    WR[3] = WR[0] * uRot[2];
                    WR[4] = p / (kapa - 1) + (WR[1] * WR[1] + WR[2] * WR[2] + WR[3] * WR[3]) / (2 * WR[0]) - WR[0] * Math.pow(Mat.L2Norm(Mat.cross(omega, elem.currentX)), 2) / 2;
                    WR[5] = 0;
                } else {
                    System.arraycopy(WL, 0, WR, 0, nEqs);
                }
                break;

            case (BoundaryType.ROTOR): // rotor
                if (isDiffusive) {
                    double[] uRot = Mat.cross(omegaRotor, Mat.minusVec(elem.currentX, rotorExcentricity));
                    WR[0] = WL[0];
                    WR[1] = WR[0] * uRot[0];
                    WR[2] = WR[0] * uRot[1];
                    WR[3] = WR[0] * uRot[2];
                    WR[4] = p / (kapa - 1) + (WR[1] * WR[1] + WR[2] * WR[2] + WR[3] * WR[3]) / (2 * WR[0]) - WR[0] * Math.pow(Mat.L2Norm(Mat.cross(omega, elem.currentX)), 2) / 2;
                    WR[5] = 0;
                } else {
                    System.arraycopy(WL, 0, WR, 0, nEqs);
                }
                break;

            case (BoundaryType.INLET):
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
                    WR[5] = Rinl * vtIn;
                } else { // supersonicky vstup
                    WR[0] = WIn[0];
                    WR[1] = WIn[1];
                    WR[2] = WIn[2];
                    WR[3] = WIn[3];
                    WR[4] = WIn[4];
                    WR[5] = WIn[0] * WIn[5];
                }
                break;

            case (BoundaryType.OUTLET):
                double ro = WL[0];
                double uo = WL[1] / WL[0];
                double vo = WL[2] / WL[0];
                double wo = WL[3] / WL[0];
                double vtOut = WL[5] / WL[0];
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
                WR[5] = ro * vtOut;
                break;

            case (BoundaryType.INVISCID_WALL): // nevazka stena
                System.arraycopy(WL, 0, WR, 0, nEqs);
                break;
        }
        return WR;
    }

    public double max(double a, double b) {
        if (a > b) {
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
        if (TT == BoundaryType.WALL || TT == BoundaryType.STATOR || TT == BoundaryType.ROTOR) {
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
                double xi = Math.max(W[dim+2],0); // vt/v
                double fv1 = xi * xi * xi / (xi * xi * xi + cv1 * cv1 * cv1);
                return new double[]{W[dim + 2]*fv1};
                
            default:
                throw new UnsupportedOperationException("undefined value " + name);
        }
    }
}
