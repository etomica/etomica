/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

//
//Ab initio potential energy surface for the nitrogen molecule pair
//
//          and thermophysical properties of nitrogen gas
//
//
//
//                           Supplement
//
//
//
//This document contains:
//
//
//
//- all details of the two N2-N2 potential functions
//
//- all 408 calculated interaction energies at different levels of theory
//
//- calculated values for the second and third pressure virial coefficients
//
//- calculated values for the zero-density limits of shear viscosity and thermal conductivity
//
//
//
//It should be viewed with a monospaced font (e.g. Courier).
//
//
//
//For questions: robert.hellmann@uni-rostock.de
//
//Look at /usr/users/rsubrama/workspace/Etomica/P2NitrogenHellmann_supplement.txt for further information
//
//---------------------------------------------------------------------------------------------------
//
//
//
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;

import etomica.space3d.Space3D;
import etomica.units.Degree;
import etomica.units.Kelvin;
import etomica.util.RandomMersenneTwister;

import etomica.util.RandomNumberGeneratorUnix;


public class P2NitrogenHellmann implements IPotentialAtomic {
    public static void main(String[] args) {
        ISpace space = Space3D.getInstance();
        int [] seeds = RandomNumberGeneratorUnix.getRandSeedArray();        
        IRandom rand1 = new RandomMersenneTwister(seeds);
        P2NitrogenHellmann pN2 = new P2NitrogenHellmann(space, rand1);        
        FileReader fileReader = null;
        String fileName = "P2NitrogenHellmann_energies.dat";
        double [] r12 = new double [408];
        double [] th1 = new double [408];
        double [] th2 = new double [408];
        double [] phi = new double [408];
        double [] eValues = new double [408];
        
        try {
            fileReader = new FileReader(fileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);            
            for (int i=0; i < 408; i++) {
                String[] str = bufReader.readLine().trim().split(" +");
                r12[i] = Double.valueOf(str[0]).doubleValue();
                th1[i] = Degree.UNIT.toSim(Double.valueOf(str[1]).doubleValue());
                th2[i] = Degree.UNIT.toSim(Double.valueOf(str[2]).doubleValue());
                phi[i] = Degree.UNIT.toSim(Double.valueOf(str[3]).doubleValue());
                eValues[i] = Double.valueOf(str[4]).doubleValue();                
            }            
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }        
        
        for (int i=0; i < 408; i++) {
            pN2.vN2Angles(r12[i], th1[i], th2[i], phi[i]);
        }
        
    }
    
    protected IBoundary boundary;
    
    protected double[][] A,alpha,b,c6;
    protected static final double[] q = {-832.77884541,1601.24507755,-1536.93246428,1601.24507755, -832.77884541};
    protected static final double[] l = {0, 0.447763006688, 0.680065710389};
    protected final IRandom random;
    protected final ISpace space;
    
    
    public P2NitrogenHellmann(ISpace space, IRandom random) {
        this.random = random;
        this.space = space;
        A = new double[5][5];
        alpha = new double[5][5];
        b = new double[5][5];
        c6 = new double[5][5];
        
        fillData();
    }
       
    public double vN2Angles (double R12, double th1, double th2, double phi) {
        /* Method to calculate the potential between site i0 on molecule 0
         * and site i1 of molecule 1 for a linear diatomic molecule with
         * 5 sites on each molecule.
         * 
         * This method takes as input, the center of mass distance and 
         * orientations of the two molecules represented using the angles
         * theta_0 , theta_1 and phi. We will almost never use this method
         * except for dubugging purposes.
         * 
         * See http://dx.doi.org/10.1080/00268976.2012.726379 for a similar 
         * definition of different angles. 
         * 
         * See http://dx.doi.org/10.1063/1.2975220 for a similar pictorial
         * representation.
         * 
         * R12  = Center of mass distance in Angstroms
         * th0 = Angle theta_0 in radians
         * th1 = Angle theta_1 in radians
         * phi = Angle phi in radians
         *   
         */

        
        IVectorMutable ex = space.makeVector();
        IVectorMutable ey = space.makeVector();
        IVectorMutable ez = space.makeVector();
        IVectorMutable a0 = space.makeVector();
        IVectorMutable a1 = space.makeVector();        
        IVectorMutable site0 = space.makeVector();
        IVectorMutable site1 = space.makeVector();
        IVectorMutable dr = space.makeVector();
                
        double cth1 = Math.cos(th1);
        double cth2 = Math.cos(th2);
        double cphi2 = Math.cos(phi);
        double sth1 = Math.sin(th1);
        double sth2 = Math.sin(th2);
        double sphi2 = Math.sin(phi);
        
        ex.E(0);
        ex.setX(0, 1);
        ey.E(0);
        ey.setX(1, 1);
        ez.E(0);
        ez.setX(2, 1);        
        
        a0.E(ex);
        a1.E(ex);
        boolean flag1 = false;
        boolean flag2 = false;
        if (th1 == 0) {
            flag1 = true;
        }
        else {
            rotateBy(cth1,sth1,ez,a0);
        }
        if (th2 == 0) {
            flag2 = true;            
        }
        else {
            rotateBy(cth2,sth2,ez,a1);
        }
        if (phi != 0) {
            if (!flag2) {
                rotateBy(cphi2,sphi2,ex,a1);
            }
            else if (!flag1) {
                rotateBy(cphi2,sphi2,ex,a0);
            }
        }
        
        double dth1 = Degree.UNIT.fromSim(Math.acos(a0.dot(ex)));
        double dth2 = Degree.UNIT.fromSim(Math.acos(a1.dot(ex)));
        IVectorMutable n0 = space.makeVector();
        IVectorMutable n1 = space.makeVector();
        n0.E(a0);
        n0.PEa1Tv1(-cth1, ex);
        if (n0.isZero())  n0.E(ey);
        n0.normalize();

        n1.E(a1);        
        n1.PEa1Tv1(-cth2, ex);
        if (n1.isZero())  n1.E(ey);
        n1.normalize();
        
        
        if (n0.isNaN() || n1.isNaN()) throw new RuntimeException("oops");

        double cphi = n0.dot(n1);
        if (cphi > 1.0) cphi = 1.0;
        if (cphi < -1.0) cphi = -1.0;
        double dphi = Degree.UNIT.fromSim(Math.acos(cphi));
        
        double v = 0;
        dr.E(0);
        dr.setX(0, R12);        
        for (int i0 = 0; i0 < 5; i0++) {
            site0.E(0);            
            if (i0 > 2) {
                site0.PEa1Tv1(l[i0 - 2], a0);
            }
            else {
                site0.PEa1Tv1(-l[2 - i0], a0);
            }            

            for (int i1 = 0; i1 < 5; i1++) {
                site1.E(dr);
                if (i1 > 2) {
                    site1.PEa1Tv1(l[i1 - 2], a1);
                }
                else {
                    site1.PEa1Tv1(-l[2 - i1], a1);
                }
                double r2ij = Math.sqrt(site0.Mv1Squared(site1));                
                double term1 = A[i0][i1]*Math.exp(-alpha[i0][i1]*r2ij);
                double r6 = r2ij*r2ij*r2ij*r2ij*r2ij*r2ij;
                double term2 = -f6(b[i0][i1]*r2ij)*c6[i0][i1]/r6;
                double term3 = q[i0]*q[i1]/r2ij;
                v += Kelvin.UNIT.toSim(term1+term2+term3);
                if (Double.isNaN(v)) throw new RuntimeException("oops "+v);
            }
        }
//        System.out.println(R12+" "+dth1 + " "+dth2+" "+dphi+" "+Kelvin.UNIT.fromSim(v));
        return v;
    }
    
    public double vN2Vectors (double R12, IVectorMutable or0, IVectorMutable or1) {
        /* Method to calculate the potential between site i0 on molecule 0
         * and site i1 of molecule 1 for a linear diatomic molecule with
         * 5 sites on each molecule.
         * 
         * This method takes as input, the center of mass distance and 
         * orientations of the two molecules represented using the vectors
         * or0 and or1  
         * 
         * R12 = Center of mass distance in Angstroms
         * or0 = unit vector pointing in the direction of orientation of 
         * molecule 0
         * or1 = unit vector pointing in the direction of orientation of
         * molecule 1
         * i0 = index of site in molecule 0
         * i1 = index of site in molecule 1  
         */
        IVectorMutable dr = space.makeVector();
        IVectorMutable a0 = space.makeVector();
        IVectorMutable a1 = space.makeVector();
        IVectorMutable site0 = space.makeVector();
        IVectorMutable site1 = space.makeVector();
                
        a0.E(or0);
        a1.E(or1);
        
        dr.E(0);
        dr.setX(0, R12);
        
        double v = 0;
        for (int i0 = 0; i0 < 5; i0++) {
            site0.E(0);
            if (i0 > 2) {
                site0.PEa1Tv1(l[i0 - 2], a0);
            }
            else {
                site0.PEa1Tv1(-l[2 - i0], a0);
            }            

            for (int i1 = 0; i1 < 5; i1++) {
                site1.E(dr);
                if (i1 > 2) {
                    site1.PEa1Tv1(l[i1 - 2], a1);
                }
                else {
                    site1.PEa1Tv1(-l[2 - i1], a1);
                }               
                double r2ij = Math.sqrt(site0.Mv1Squared(site1));                
                double term1 = A[i0][i1]*Math.exp(-alpha[i0][i1]*r2ij);
                double r6 = r2ij*r2ij*r2ij*r2ij*r2ij*r2ij;
                double term2 = -f6(b[i0][i1]*r2ij)*c6[i0][i1]/r6;
                double term3 = q[i0]*q[i1]/r2ij;
                v += Kelvin.UNIT.toSim(term1+term2+term3);
                if (Double.isNaN(v)) throw new RuntimeException("oops "+v);
            }
        }        
        return v;
    }
    
    protected double f6 (double bR) {
        double term = 1;
        double sum = 1;
        for (int i = 1; i <= 6; i++) {
            term *= bR/(double)i;
            sum += term;
        }
        sum *= -Math.exp(-bR);
        sum += 1.0;
        if (Double.isNaN(sum) || Double.isNaN(term) || Double.isInfinite(sum) || Double.isInfinite(term)) throw new RuntimeException(" oops!"+sum+term);
        return sum;
    }

    protected void fillData() {                
        FileReader fileReader = null;
        String fileName = "P2NitrogenHellmann_Aij.dat";
        
        try {
            fileReader = new FileReader(fileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);            
            for (int i=0; i<5; i++) {
                String[] str = bufReader.readLine().trim().split(" +");
                for (int j=0; j<5; j++) {                    
                    A[i][j] = Double.valueOf(str[j]).doubleValue();
                }
            }

            
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        
        fileName = "P2NitrogenHellmann_alphaij.dat";
        
        try {
            fileReader = new FileReader(fileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);            
            for (int i=0; i<5; i++) {
                String[] str = bufReader.readLine().trim().split(" +");
                for (int j=0; j<5; j++) {                    
                    alpha[i][j] = Double.valueOf(str[j]).doubleValue();
                }
            }

            
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        
        fileName = "P2NitrogenHellmann_bij.dat";
        
        try {
            fileReader = new FileReader(fileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);            
            for (int i=0; i<5; i++) {
                String[] str = bufReader.readLine().trim().split(" +");
                for (int j=0; j<5; j++) {                    
                    b[i][j] = Double.valueOf(str[j]).doubleValue();
                }
            }

            
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        
        fileName = "P2NitrogenHellmann_c6ij.dat";
        
        try {
            fileReader = new FileReader(fileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);            
            for (int i=0; i<5; i++) {
                String[] str = bufReader.readLine().trim().split(" +");
                for (int j=0; j<5; j++) {                    
                    c6[i][j] = Double.valueOf(str[j]).doubleValue();
                }
            }            
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
    }
    public static void rotateBy(double cdt, double sdt, IVector axis, IVectorMutable direction) {
        // consider a circle on the surface of the unit sphere.  The given axis
        // passes through the center of the circle.  The circle passes through
        // the current direction vector and the vector v4 defined below.  We
        // rotate the direction by the given angle (dt) around the circle.
        
        // v1 is the projection of direction onto axis
        // v2 is the component of direction perpendicular to axis
        // v3 has the same magnitude as v2 and is perpendicular to both
        //    direction and axis
        // v4 is a unit vector whose components are v1 and v3
        
        // v1 = v1overAxis * axis
        ISpace space = Space3D.getInstance();        
        double v1overAxis = axis.dot(direction);
        IVectorMutable temp = space.makeVector();
        IVectorMutable temp2 = space.makeVector();
        temp.Ea1Tv1(-v1overAxis, axis);
        temp.PE(direction);
        // now temp = v2
        temp2.E(axis);
        temp2.XE(direction);
        // now temp2 = v3
        direction.Ea1Tv1(cdt, temp);
        direction.PEa1Tv1(sdt, temp2);
        direction.PEa1Tv1(v1overAxis, axis);
        
    }
    public double getRange() {        
        return Double.POSITIVE_INFINITY;
    }
    
    public void setBox(IBox box) {
        boundary = box.getBoundary();    
    }

    public int nBody() {
        return 2;
    }   
     
    public double energy(IAtomList atoms) {
    
        return 0;
    }
}

