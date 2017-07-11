/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import java.io.FileWriter;
import java.io.IOException;

import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;

public class P2HydrogenPatkowskiIso extends Potential2SoftSpherical {

    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        P2HydrogenPatkowskiIso pot1 = new P2HydrogenPatkowskiIso(space);
        double r2limit = BohrRadius.UNIT.toSim(4.75)*BohrRadius.UNIT.toSim(4.75);
        try{                        
            FileWriter one = new FileWriter("core.dat");           
            int points = 1000;
            for (int i = 1; i<points; i++) {
                double r2 = i*r2limit/points;
                double u1 = pot1.u(r2);
                double du1 = pot1.du(r2)/(Math.sqrt(r2));
                double d2u1 = pot1.d2u(r2)/r2;
                one.write(Math.sqrt(r2)+" "+u1+" "+du1+" "+d2u1+"\n");
//                one.write(Math.sqrt(r2)+" "+du1+" "+d2u1+"\n");
            }
            one.close();
        }        
        catch(IOException ex1){
            throw new RuntimeException(ex1);
        }
        System.out.println("Energy = "+pot1.u(r2limit));
    }

    public P2HydrogenPatkowskiIso(Space space) {
        super(space);                
        // TODO Auto-generated constructor stub
    }

    public double u(double r2) {        
        //conversion factor (K -> a.u.) (1 K = 3.16669^(-6) a.u.)
        double xK2au = 3.16669E-6;
        double Econv= 1.0/xK2au;
        double R = BohrRadius.UNIT.fromSim(Math.sqrt(r2));
        if (R < 1.1) return Double.POSITIVE_INFINITY; // setting the core at 1.1 bohrs
        r2 = R*R;
        //compute the long range terms
        double rr02 = 1.0/r2;        
        double rr06 = rr02*rr02*rr02;
        double rr08 = rr06*rr02;
        double rr10 = rr08*rr02;
        double elong06 = cdataIso[0]*rr06*Econv;
        double elong08 = cdataIso[1]*rr08*Econv;
        double elong10 = cdataIso[2]*rr10*Econv;

        //use the damping factors
        double damp = -2.0*cexIso[1];
        double ddd06 = d(6,damp,R);
        double ddd08 = d(8,damp,R);
        double ddd10 = d(10,damp,R);
        double vlong = ddd06*elong06 + ddd08*elong08 + ddd10*elong10;        

        //compute the short-range part
        double vex = Math.exp(2.0*(cexIso[0] + cexIso[1]*R));        
        double vsp = 2.0*(cspIso[0] + cspIso[1]*R + cspIso[2]*R*R + cspIso[3]*R*R*R);
        double vshort = vex*vsp;        

        //sum up        
        double en = vshort+vlong;   
        return en;
    }
    protected double d(int n,double beta, double r) {
        //
        //             calculate the damping factor (small R correct)
        //
        
        double br=beta*r;
        double sum=1.00;
        double term=1.00;          
        for (int i=1; i<=n; i++) {
            term=term*br/i;
            sum=sum+term;
        }
        double d1 = 1.00 - Math.exp(-br)*sum;
        if(Math.abs(d1) < 1.0e-8) {                        
            d1=0.00;
            for (int i=n+1; i<=1000; i++) {
                term=term*br/i;
                d1=d1+term;
                if(term/d1 < 1.0e-8) break;
            }
            if(term/d1 >= 1.0e-8) throw new RuntimeException("No convergence in d");
            d1=d1*Math.exp(-br);
        }        
        return d1;
    }    

    public double du(double r2) {

        double xK2au = 3.16669E-6;
        double Econv= 1.0/xK2au;
        double R = BohrRadius.UNIT.fromSim(Math.sqrt(r2));
        r2 = R*R;
        double rr02 = 1.0/r2;        
        double rr06 = rr02*rr02*rr02;
        double rr08 = rr06*rr02;
        double rr10 = rr08*rr02;

        double damp = -2.0*cexIso[1];
        double ddd06 = d(6,damp,R);
        double ddd08 = d(8,damp,R);
        double ddd10 = d(10,damp,R);

        double elong06 = cdataIso[0]*rr06*Econv;
        double elong08 = cdataIso[1]*rr08*Econv;
        double elong10 = cdataIso[2]*rr10*Econv;
        double delong06dr = -6*elong06/R;
        double delong08dr = -8*elong08/R;
        double delong10dr = -10*elong10/R;

        double vex = Math.exp(2.0*(cexIso[0] + cexIso[1]*R));        
        double dvexdr = 2.0*cexIso[1]*vex;        
        double vsp = 2.0*(cspIso[0] + cspIso[1]*R + cspIso[2]*R*R + cspIso[3]*R*R*R);        
        double dvspdr = 2.0*(cspIso[1] + 2*cspIso[2]*R + 3*cspIso[3]*R*R);
        double dvshortdr = vex*dvspdr + vsp*dvexdr;
//        dvshortdr = 0.0;
        

        double d06dr = 0.0;
        double d08dr = 0.0;
        double d10dr = 0.0;
        double term0 = 1.00;
        double sum = 1.00;

        for (int i = 1; i <= 6;i++) {
            term0 *= damp*R/i;// (damp*R)^6/6!
            sum += term0;
        }
        double d1 = 1.00 - Math.exp(-damp*R)*sum;
        if (Math.abs(d1) < 1.0E-8){
            d1 = 0.00;
            double term1 = term0;
            for (int i = 7; i <= 1000; i++) {
                term1 *= damp*R/i;
                d1 += term1;
                if (term1/d1 < 1.0E-8) break;
            }
            if (term1/d1 >= 1.0E-8) throw new RuntimeException("No convergence in d");           
            d06dr = damp*Math.exp(-damp*R)*(term0 - term1);
        }
        else {
            d06dr = term0*damp*Math.exp(-damp*R);
        }
        
        term0 = 1.00;
        sum = 1.00;
        for (int i = 1; i <= 8;i++) {
            term0 *= damp*R/i;
            sum += term0;
        }
        d1 = 1.00 - Math.exp(-damp*R)*sum;
        if (Math.abs(d1) < 1.0E-8) {            
            d1 = 0.00;
            double term1 = term0;
            for (int i = 9; i <= 1000; i++) {
                term1 *= damp*R/i;
                d1 += term1;
                if (term1/d1 < 1.0E-8) break;
            }
            if (term1/d1 >= 1.0E-8) throw new RuntimeException("No convergence in d");
            d08dr = damp*Math.exp(-damp*R)*(term0 - term1);
        }
        else {
            d08dr = term0*damp*Math.exp(-damp*R);
        }

        term0 = 1.00;
        sum = 1.00;
        for (int i = 1; i <= 10;i++) {
            term0 *= damp*R/i;
            sum += term0;
        }
        d1 = 1.00 - Math.exp(-damp*R)*sum;
        if (Math.abs(d1) < 1.0E-8) {            
            d1 = 0.00;
            double term1 = term0;
            for (int i = 11; i <= 1000; i++) {
                term1 *= damp*R/i;
                d1 += term1;
                if (term1/d1 < 1.0E-8) break;
            }
            if (term1/d1 >= 1.0E-8) throw new RuntimeException("No convergence in d");
            d10dr = damp*Math.exp(-damp*R)*(term0 - term1);
        }
        else {
            d10dr = term0*damp*Math.exp(-damp*R);
        }

        double dvlongdr = ddd06*delong06dr + elong06*d06dr + ddd08*delong08dr + elong08*d08dr + ddd10*delong10dr+ elong10*d10dr;        
        double rDudr = R*(dvshortdr + dvlongdr);
        return rDudr;
    }

    public double d2u(double r2) {

        double xK2au = 3.16669E-6;
        double Econv= 1.0/xK2au;
        double R = BohrRadius.UNIT.fromSim(Math.sqrt(r2));
        r2 = R*R;
        double rr02 = 1.0/r2;        
        double rr06 = rr02*rr02*rr02;
        double rr08 = rr06*rr02;
        double rr10 = rr08*rr02;
        
        double vex = Math.exp(2.0*(cexIso[0] + cexIso[1]*R));        
        double vsp = 2.0*(cspIso[0] + cspIso[1]*R + cspIso[2]*R*R + cspIso[3]*R*R*R);        
        double dvexdr = 2.0*cexIso[1]*vex;        
        double dvspdr = 2.0*(cspIso[1] + 2*cspIso[2]*R + 3*cspIso[3]*R*R);        
        double d2vexdr2 = 2.0*cexIso[1]*dvexdr;        
        double d2vspdr2 = 2.0*(2*cspIso[2]+6*cspIso[3]*R);
        double d2vshortdr2 = vex*d2vspdr2 + 2*dvexdr*dvspdr + vsp*d2vexdr2;
//        d2vshortdr2 = 0.0;

        double damp = -2.0*cexIso[1];
        double ddd06 = d(6,damp,R);
        double ddd08 = d(8,damp,R);
        double ddd10 = d(10,damp,R);
        double d06dr = 0.0;
        double d08dr = 0.0;
        double d10dr = 0.0;
        double d206dr2 = 0.0;
        double d208dr2 = 0.0;
        double d210dr2 = 0.0;
        double term0 = 1.00;
        double sum = 1.00;
        
        for (int i = 1; i <= 6;i++) {
            term0 *= damp*R/i;// (damp*R)^6/6!
            sum += term0;
        }
        double d1 = 1.00 - Math.exp(-damp*R)*sum;
        if (Math.abs(d1) < 1.0E-8){
            int m = 0;
            d1 = 0.00;
            double term1 = term0;
            for (int i = 7; i <= 1000; i++) {
                term1 *= damp*R/i;
                d1 += term1;
                if (term1/d1 < 1.0E-8) {
                    m = i;
                    break;
                }
            }
            if (term1/d1 >= 1.0E-8) throw new RuntimeException("No convergence in d");           
            d06dr = damp*Math.exp(-damp*R)*(term0 - term1);
            d206dr2 = damp*Math.exp(-damp*R)*(term0*6/R - term1*m/R) - damp*d06dr;
        }
        else {
            d06dr = term0*damp*Math.exp(-damp*R);
            d206dr2 = d06dr*(6/R - damp);
        }
        
        term0 = 1.00;
        sum = 1.00;
        for (int i = 1; i <= 8;i++) {
            term0 *= damp*R/i;
            sum += term0;
        }
        d1 = 1.00 - Math.exp(-damp*R)*sum;
        if (Math.abs(d1) < 1.0E-8) {
            int m = 0;
            d1 = 0.00;
            double term1 = term0;
            for (int i = 9; i <= 1000; i++) {
                term1 *= damp*R/i;
                d1 += term1;
                if (term1/d1 < 1.0E-8) {
                    m = i;
                    break;
                }
            }
            if (term1/d1 >= 1.0E-8) throw new RuntimeException("No convergence in d");
            d08dr = damp*Math.exp(-damp*R)*(term0 - term1);
            d208dr2 = damp*Math.exp(-damp*R)*(term0*8/R - term1*m/R) - damp*d08dr;
        }
        else {
            d08dr = term0*damp*Math.exp(-damp*R);
            d208dr2 = d08dr*(8/R - damp);
        }

        term0 = 1.00;
        sum = 1.00;
        for (int i = 1; i <= 10;i++) {
            term0 *= damp*R/i;
            sum += term0;
        }
        d1 = 1.00 - Math.exp(-damp*R)*sum;
        if (Math.abs(d1) < 1.0E-8) {            
            int m = 0;
            d1 = 0.00;
            double term1 = term0;
            for (int i = 11; i <= 1000; i++) {
                term1 *= damp*R/i;
                d1 += term1;
                if (term1/d1 < 1.0E-8) {
                    m = i;
                    break;
                }
            }
            if (term1/d1 >= 1.0E-8) throw new RuntimeException("No convergence in d");
            d10dr = damp*Math.exp(-damp*R)*(term0 - term1);
            d210dr2 = damp*Math.exp(-damp*R)*(term0*10/R - term1*m/R) - damp*d10dr;
        }
        else {
            d10dr = term0*damp*Math.exp(-damp*R);
            d210dr2 = d10dr*(10/R - damp);
        }
        
        double elong06 = cdataIso[0]*rr06*Econv;
        double elong08 = cdataIso[1]*rr08*Econv;
        double elong10 = cdataIso[2]*rr10*Econv;
        double delong06dr = -6*elong06/R;
        double delong08dr = -8*elong08/R;
        double delong10dr = -10*elong10/R;
        double d2elong06dr2 = -7*delong06dr/R;
        double d2elong08dr2 = -9*delong08dr/R;
        double d2elong10dr2 = -11*delong10dr/R;        
//        double dvlongdr = ddd06*delong06dr + elong06*d06dr + ddd08*delong08dr + elong08*d08dr + ddd10*delong10dr+ elong10*d10dr;

        double d2vlongdr2 = ddd06*d2elong06dr2 + 2*delong06dr*d06dr + elong06*d206dr2 + ddd08*d2elong08dr2 + 2*delong08dr*d08dr + elong08*d208dr2 + ddd10*d2elong10dr2 + 2*d10dr*delong10dr + elong10*d210dr2;       
        double r2d2udr2 = r2*(d2vshortdr2 + d2vlongdr2);
        return r2d2udr2;
    }

    public double uInt(double rC) {
        // TODO Auto-generated method stub
        return 0;
    }
    protected static final double [] cexIso = new double [] {0.79181110538811E+01,-0.86728533815761E+00};
    protected static final double [] cspIso = new double [] {-0.13759109019103E+00, 0.18174129593930E+00,-0.27161233180350E-01,0.10384349156326E-02};
    protected static final double [] cdataIso = new double [] {-12.058168,-213.6,-4700.0};    
}

