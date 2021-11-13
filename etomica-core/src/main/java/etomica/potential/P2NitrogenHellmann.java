/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.chem.elements.Nitrogen;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.Degree;
import etomica.units.Kelvin;
import etomica.util.Constants;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

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


public class P2NitrogenHellmann implements Potential2Soft {
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        P2NitrogenHellmann pN2 = new P2NitrogenHellmann(space);
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
    
    protected final double[][] A,alpha,b,c6;
    protected static final int [] siteID = {0,1,2,1,0};
    protected double[] q, pos;
    public boolean parametersB = true; // set this to false if using parameters for V_12^A potential
    protected static final double[] AA = {0.846144529794E7, -0.825811757649E7, 0.116781452519E8, 0.198375106632E7, -0.889934931566E7, 0.216444280602E8};
    protected static final double[] alphaA = {2.93660213026, 2.64088584245, 2.92688055782, 2.11392143479, 3.05000862940, 3.22097240942};
    protected static final double[] bA = {3.10077160627, 3.52966064535, 2.92682684017, 2.44007490969, 1.97080180486, 1.89759495054};
    protected static final double[] qA = {-794.51215612, 1621.86520764, -1654.70610304, 1621.86520764, -794.51215612};
    protected static final double[] c6A = {0.224212900521E7, -0.438154598247E7, 0.304119164146E7, 0.103132819062E8, -0.740002575021E7, 0.277494715937E7};
    protected static final double[] sitePosA = {-0.682390307412, -0.434260425799, 0.00, 0.434260425799, 0.682390307412};
    protected static final double[] AB = {0.973347918383E7,-0.954555977809E7,0.122158259267E8,0.299460243665E7,-0.819908034347E7, 0.163947777734E8};
    protected static final double[] alphaB = {3.06144571072,2.58710992361,2.96686681629,2.15319940621,2.84661195657,2.99548316813};
    protected static final double[] bB = {2.58031350518, 3.45760438302, 2.46746232590, 2.42577961527, 2.02508542307, 1.97117981681};
    protected static final double[] qB = {-832.77884541,1601.24507755,-1536.93246428,1601.24507755, -832.77884541};
    protected static final double[] c6B = {0.298807116692E7, -0.608284467163E7, 0.490318811890E7, 0.146889670654E8, -0.129841807274E8, 0.107874613877E8};
    protected static final double[] sitePosB = {-0.680065710389,-0.447763006688, 0.00, 0.447763006688, 0.680065710389};
    protected static final double dHSCore = 2.0;
    protected final Space space;
    protected static final double massN2 = 2*Nitrogen.INSTANCE.getMass();
    public static final double blN2 = 1.1014;
    protected static final double moment = 0.25*massN2*blN2*blN2;
    
    public P2NitrogenHellmann(Space space) {
        this.space = space;
        A = new double[3][3];
        alpha = new double[3][3];
        b = new double[3][3];
        c6 = new double[3][3];
        q = new double[3];
        pos = new double[5];    
        gradientAndTorque = new Vector[2][2];
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                gradientAndTorque[i][j] = space.makeVector();
            }
        }
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

        if (R12 < dHSCore) return Double.POSITIVE_INFINITY;
        
        Vector ex = space.makeVector();
        Vector ey = space.makeVector();
        Vector ez = space.makeVector();
        Vector a0 = space.makeVector();
        Vector a1 = space.makeVector();
        Vector site0 = space.makeVector();
        Vector site1 = space.makeVector();
        Vector dr = space.makeVector();
                
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
        Vector n0 = space.makeVector();
        Vector n1 = space.makeVector();
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
        for (int i = 0; i < 5; i++) {
            int i0 = siteID[i];
            site0.E(0);
            site0.PEa1Tv1(pos[i], a0);
            for (int j = 0; j < 5; j++) {
                int i1 = siteID[j];
                site1.E(dr);
                site1.PEa1Tv1(pos[j], a1);                
                double rij = Math.sqrt(site0.Mv1Squared(site1));                                
                double term1 = A[i0][i1]*Math.exp(-alpha[i0][i1]*rij);
                double r6 = rij*rij*rij*rij*rij*rij;
                double term2 = -f6(b[i0][i1]*rij)*c6[i0][i1]/r6;
                double term3 = q[i0]*q[i1]/rij;
                v += term1 + term2 + term3;
                if (Double.isNaN(v)) throw new RuntimeException("oops "+v);
            }
        }
//        System.out.println(R12+" "+dth1 + " "+dth2+" "+dphi+" "+v);
        return Kelvin.UNIT.toSim(v);
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
        int k = 0;
        if (parametersB) {
            for (int i=0; i < 3; i++) {
                A[i][i] = AB[k];
                alpha[i][i] = alphaB[k];
                b[i][i] = bB[k];
                c6[i][i] = c6B[k];            
                k++;
                for (int j=i+1; j < 3; j++) {
                    A[i][j] = A[j][i] = AB[k];
                    alpha[i][j] = alpha[j][i] = alphaB[k];
                    b[i][j] = b[j][i] = bB[k];
                    c6[i][j] = c6[j][i] = c6B[k];                
                    k++;
                }
            }
            for (int i = 0; i < q.length; i++) {
                q[i] = qB[i];
            }
            pos = sitePosB;
        }
        else {
            for (int i=0; i < 3; i++) {
                A[i][i] = AA[k];
                alpha[i][i] = alphaA[k];
                b[i][i] = bA[k];
                c6[i][i] = c6A[k];            
                k++;
                for (int j=i+1; j < 3; j++) {
                    A[i][j] = A[j][i] = AA[k];
                    alpha[i][j] = alpha[j][i] = alphaA[k];
                    b[i][j] = b[j][i] = bA[k];
                    c6[i][j] = c6[j][i] = c6A[k];                
                    k++;
                }
            }
            for (int i = 0; i < q.length; i++) {
                q[i] = qA[i];
            }
            pos = sitePosA;
        }
    }
    
    public static void rotateBy(double cdt, double sdt, Vector axis, Vector direction) {
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
        Space space = Space3D.getInstance();
        double v1overAxis = axis.dot(direction);
        Vector temp = space.makeVector();
        Vector temp2 = space.makeVector();
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

    public int nBody() {
        return 2;
    }   
     
    public double energy(IAtomList atoms) {
        IAtomOriented a0 = (IAtomOriented)atoms.get(0);
//        int aIndex0 = a0.getLeafIndex();
        IAtomOriented a1 = (IAtomOriented)atoms.get(1);
//        int aIndex1 = a1.getLeafIndex();
        Vector or0 = space.makeVector();
        or0.E(a0.getOrientation().getDirection());
        Vector or1 = space.makeVector();
        or1.E(a1.getOrientation().getDirection());                
        Vector com0 = a0.getPosition();
        Vector com1 = a1.getPosition();
//        double R12 = Math.sqrt(com0.Mv1Squared(com1));
//        if (aIndex0 == 1 && aIndex1 == 2) {
//            IVectorMutable p1 = space.makeVector();
//            p1.E(com0);
//            p1.PEa1Tv1(pos[4],hh0);
//            double r12 = Math.sqrt(com1.Mv1Squared(p1));
//            if (r12 < dHSCore) return Double.POSITIVE_INFINITY;
//            return 0;
//        }
        Vector ex = space.makeVector();
        Vector ey = space.makeVector();
        Vector ez = space.makeVector();
                
        Vector site0 = space.makeVector();
        Vector site1 = space.makeVector();
        
        ex.E(0);
        ex.setX(0, 1);
        ey.E(0);
        ey.setX(1, 1);
        ez.E(0);
        ez.setX(2, 1);
        
        double th1 = Math.acos(or0.dot(ex));
        double th2 = Math.acos(-or1.dot(ex));
        double dth1 = Degree.UNIT.fromSim(th1);
        double dth2 = Degree.UNIT.fromSim(th2);
        double cth1 = Math.cos(Degree.UNIT.toSim(dth1));
        double cth2 = Math.cos(Degree.UNIT.toSim(dth2));
        Vector normal0 = space.makeVector();
        Vector normal1 = space.makeVector();
        normal0.E(or0);
        normal0.PEa1Tv1(-cth1, ex);
        if (normal0.isZero())  normal0.E(ey);
        normal0.normalize();

        normal1.E(or1);        
        normal1.PEa1Tv1(cth2, ex);
        if (normal1.isZero()) normal1.E(ey);
        normal1.normalize();
        
        if (normal0.isNaN() || normal1.isNaN()) throw new RuntimeException("oops");
        double cphi = normal0.dot(normal1);
        if (cphi > 1.0) cphi = 1.0;
        if (cphi < -1.0) cphi = -1.0;
        double phi = Math.acos(cphi);
        if (th1 > Math.PI/2.0) {
            th1 = Math.PI - th1;            
        }
        if (th2 > Math.PI/2.0) {
            th2 = Math.PI - th2;            
        }
        if (phi > Math.PI) {
            phi = 2*Math.PI - phi;
        }
        double dphi = Degree.UNIT.fromSim(phi);
        double v = 0;
        double R12 = Math.sqrt(com0.Mv1Squared(com1));
        if (R12 < dHSCore) return Double.POSITIVE_INFINITY;
        for (int i = 0; i < 5; i++) {
            int i0 = siteID[i];
            site0.E(com0);
            site0.PEa1Tv1(pos[i], or0);
            for (int j = 0; j < 5; j++) {
                int i1 = siteID[j];
                site1.E(com1);
                site1.PEa1Tv1(pos[j], or1);                               
                double rij = Math.sqrt(site0.Mv1Squared(site1));                
                double term1 = A[i0][i1]*Math.exp(-alpha[i0][i1]*rij);
                double r6 = rij*rij*rij*rij*rij*rij;
                double term2 = -f6(b[i0][i1]*rij)*c6[i0][i1]/r6;
                double term3 = q[i0]*q[i1]/rij;
                v += term1 + term2 + term3;
                if (Double.isNaN(v)) throw new RuntimeException("oops "+v);
            }
        }
//        System.out.println(R12+" "+dth1 + " "+dth2+" "+dphi+" "+v);        
        double E = Kelvin.UNIT.toSim(v);
        return E;
    }
    
    
    public P2N2QFH makeQFH(double temperature) {
        return new P2N2QFH(temperature);
    }
    
    public class P2N2QFH implements IPotentialAtomic {        
        protected final Vector[][] gi;
        protected final Tensor tt0Tensor, tt1Tensor, rr0Tensor, rr1Tensor;
        protected final Tensor ijTensor, rTensor0, rTensor1, identity;
        protected final Tensor ijRTensor;
        protected final Tensor rot0, rot1;
        protected final Vector or01, or11, or02, or12;
        protected final Vector[] allOr0, allOr1;
        protected final Vector drijRot;
        public double[][] d2tot = new double[2][6];
        protected final double temperature, fac;        
        
        public P2N2QFH(double temperature) { // copied from Andrew's code: P2CO2Hellmann.java
            ijTensor = space.makeTensor();
            identity = space.makeTensor();
            tt0Tensor = space.makeTensor();
            tt1Tensor = space.makeTensor();
            rr0Tensor = space.makeTensor();
            rr1Tensor = space.makeTensor();
            rTensor0 = space.makeTensor();
            rTensor1 = space.makeTensor();
            ijRTensor = space.makeTensor();
            identity.E(new double[][]{{1,0,0},{0,1,0},{0,0,1}});
            gi = new Vector[2][5];
            for (int i=0; i<5; i++) {
                gi[0][i] = space.makeVector();
                gi[1][i] = space.makeVector();
            }
            or01 = space.makeVector();
            or11 = space.makeVector();
            or02 = space.makeVector();
            or12 = space.makeVector();
            allOr0 = new Vector[]{null, or01, or02};
            allOr1 = new Vector[]{null, or11, or12};
            drijRot = space.makeVector();
            rot0 = space.makeTensor();
            rot1 = space.makeTensor();            
            
            this.temperature = temperature;
            double hbar = Constants.PLANCK_H/(2*Math.PI);
            fac = hbar*hbar/(12.00*temperature);
        }

        public double getRange() {
            return Double.POSITIVE_INFINITY;
        }

        public int nBody() {
            return 2;
        }
        protected void getPerp(Vector or, Vector perp1, Vector perp2) {
            int max = 0;
            if (Math.abs(or.getX(1)) > Math.abs(or.getX(0))) max=1;
            if (Math.abs(or.getX(2)) > Math.abs(or.getX(max))) max=2;
            int min = 0;
            if (Math.abs(or.getX(1)) < Math.abs(or.getX(0))) min=1;
            if (Math.abs(or.getX(2)) < Math.abs(or.getX(min))) min=2;
            perp1.E(0);
            perp1.setX(min, or.getX(max));
            perp1.setX(max, -or.getX(min));
            perp1.normalize();
            perp1.PEa1Tv1(-perp1.dot(or), or);
            perp1.normalize();
            perp2.E(or);
            perp2.XE(perp1);
        }

        public double energy(IAtomList atoms) {
            IAtomOriented atom0 = (IAtomOriented)atoms.get(0);
            IAtomOriented atom1 = (IAtomOriented)atoms.get(1);
            Vector cm0 = atom0.getPosition();
            Vector cm1 = atom1.getPosition();
            double R12 = Math.sqrt(cm0.Mv1Squared(cm1));
            if (R12 < dHSCore) return Double.POSITIVE_INFINITY;
            Vector or0 = atom0.getOrientation().getDirection();
            Vector site0 = space.makeVector();
            Vector site1 = space.makeVector();
            Vector dr = space.makeVector();
            getPerp(or0, or01, or02);
            allOr0[0] = or0;
            rot0.E(allOr0);
            rot0.invert();
            
            Vector or1 = atom1.getOrientation().getDirection();
            getPerp(or1, or11, or12);
            allOr1[0] = or1;
            rot1.E(allOr1);
            rot1.invert();

            tt0Tensor.E(0);
            tt1Tensor.E(0);
            rr0Tensor.E(0);
            rr1Tensor.E(0);
            for (int i=0; i<5; i++) {
                gi[0][i].E(0);
                gi[1][i].E(0);
            }
            for (int i=0; i<6; i++ ){
                d2tot[0][i] = 0;
                d2tot[1][i] = 0;
            }
            double u = 0;
            for (int i=0; i<5; i++) {
                int ii = siteID[i];
                site0.Ea1Tv1(pos[i], or0);
//                rTensor0.setComponent(0, 1, 0);
//                rTensor0.setComponent(1, 0, -0);
//                rTensor0.setComponent(0, 2, -0);
//                rTensor0.setComponent(2, 0, 0);
                rTensor0.setComponent(1, 2, pos[i]);
                rTensor0.setComponent(2, 1, -pos[i]);
                site0.PE(cm0);
                for (int j=0; j<5; j++) {
                    int jj = siteID[j];
                    site1.Ea1Tv1(pos[j], or1);
//                    rTensor1.setComponent(0, 1, 0);
//                    rTensor1.setComponent(1, 0, -0);
//                    rTensor1.setComponent(0, 2, -0);
//                    rTensor1.setComponent(2, 0, 0);
                    rTensor1.setComponent(1, 2, pos[j]);
                    rTensor1.setComponent(2, 1, -pos[j]);
                    site1.PE(cm1);
                    dr.Ev1Mv2(site1, site0);
                    double rij2 = dr.squared();
                    double rij = Math.sqrt(rij2);                    
                    double ar = alpha[ii][jj]*rij;
                    double uExp = A[ii][jj]*Math.exp(-ar);
                    double rduExpdr = -A[ii][jj]*ar*Math.exp(-ar);
                    double r2du2Expdr2 = A[ii][jj]*ar*ar*Math.exp(-ar);

                    double sum = 1;
                    double dsum = 0;
                    double d2sum = 0;
                    double br = b[ii][jj]*rij;
                    double term = 1;
                    double dterm = 1;
                    double d2term = 0;
                    for (int k=1; k<=6; k++) {
                        term *= br/k;
                        if (k==1) {
                            dterm = br;
                        }
                        else {
                            dterm *= br/(k-1);
                        }
                        if (k==2) {
                            d2term = br*br;
                        }
                        else if (k>2) {
                            d2term *= br/(k-2);
                        }
                        sum += term;
                        dsum += dterm;
                        d2sum += d2term;
                    }
                    if (sum==1) {
                        return Double.POSITIVE_INFINITY;
                    }
                    double rij6 = rij2*rij2*rij2;
                    double expbr = Math.exp(-br);
                    double u6 = -(1-expbr*sum)*c6[ii][jj]/rij6;
                    double rdu6dr = (-expbr*br*sum + expbr*dsum + 6*(1-expbr*sum))*c6[ii][jj]/rij6;
                    double r2du26dr2 = (expbr*br*br*sum-expbr*br*dsum -br*expbr*dsum+expbr*d2sum + 6*(-1 + br*expbr*sum-expbr*dsum+expbr*sum) + -6*(-expbr*br*sum + expbr*dsum + 6*(1-expbr*sum)))*c6[ii][jj]/rij6;
                    
                    double uCharge = q[ii]*q[jj]/rij;
                    double rduChargedr = -q[ii]*q[jj]/rij;
                    double r2d2uChargedr2 = 2*q[ii]*q[jj]/rij;
                    
                    u += Kelvin.UNIT.toSim(uExp + u6 + uCharge);
                    double rdudr = Kelvin.UNIT.toSim(rduExpdr + rdu6dr + rduChargedr);
                    double r2d2udr2 = Kelvin.UNIT.toSim(r2du2Expdr2 + r2du26dr2 + r2d2uChargedr2);

                    // molecule 0
                    drijRot.E(dr);
                    rot0.transform(drijRot);
                    ijTensor.Ev1v2(drijRot, drijRot);
                    ijTensor.TE((rdudr - r2d2udr2)/(rij2*rij2));
                    ijTensor.PEa1Tt1(-rdudr/rij2, identity);

                    tt0Tensor.ME(ijTensor);

                    rTensor0.transpose();
                    ijRTensor.E(rTensor0);
                    ijRTensor.TE(ijTensor);
                    rTensor0.transpose();
                    ijRTensor.TE(rTensor0);
                    rr0Tensor.ME(ijRTensor);
                    
                    drijRot.TE(rdudr/rij2);
                    gi[0][i].ME(drijRot);


                    // molecule 1
                    drijRot.E(dr);
                    rot1.transform(drijRot);
                    ijTensor.Ev1v2(drijRot, drijRot);
//                    System.out.println("r2 "+rij2+" "+ijTensor.component(0,0));
                    ijTensor.TE((rdudr - r2d2udr2)/(rij2*rij2));
                    ijTensor.PEa1Tt1(-rdudr/rij2, identity);

                    // we only need tt for 0, both have the same contribution
                    tt1Tensor.ME(ijTensor);
//                    System.out.println("d2 "+r2d2udr2/rij2+" "+ijTensor.component(0,0)+" "+tt1Tensor.component(0,0));

                    rTensor1.transpose();
                    ijRTensor.E(rTensor1);
                    ijRTensor.TE(ijTensor);
                    rTensor1.transpose();
                    ijRTensor.TE(rTensor1);
                    rr1Tensor.ME(ijRTensor);
                    
                    drijRot.TE(rdudr/rij2);
                    gi[1][j].PE(drijRot);
                }
            }
            
            for (int i=0; i<5; i++) {
                site0.E(0);
                site0.setX(0, pos[i]);
                rr0Tensor.PEv1v2(site0, gi[0][i]);
                rr0Tensor.PEa1Tt1(-pos[i]*gi[0][i].getX(0), identity);

                site1.E(0);
                site1.setX(0, pos[i]);
                rr1Tensor.PEv1v2(site1, gi[1][i]);
                rr1Tensor.PEa1Tt1(-pos[i]*gi[1][i].getX(0), identity);

            }
            double sum = 0;
            for (int i=0; i<3; i++){
                d2tot[0][i] += tt0Tensor.component(i,i);
                d2tot[1][i] += tt1Tensor.component(i,i);
                sum += tt0Tensor.component(i,i)/massN2;
            }
            for (int i=1; i<3; i++){
                d2tot[0][3+i] += rr0Tensor.component(i,i);
                d2tot[1][3+i] += rr1Tensor.component(i,i);
                sum += (rr0Tensor.component(i,i) + rr1Tensor.component(i,i))/(2*moment);
            }            
            return (u + fac*sum);
        }
    }
    protected final Vector[][] gradientAndTorque;
    
    public Vector[][] gradientAndTorque(IAtomList atoms) {
        IAtomOriented atom0 = (IAtomOriented)atoms.get(0);
        IAtomOriented atom1 = (IAtomOriented)atoms.get(1);
        Vector cm0 = atom0.getPosition();
        Vector cm1 = atom1.getPosition();
        Vector or0 = atom0.getOrientation().getDirection();
        Vector or1 = atom1.getOrientation().getDirection();
        double R12 = Math.sqrt(cm0.Mv1Squared(cm1));        
        if (R12 < dHSCore) {
            gradientAndTorque[0][0].E(Double.NaN);
            gradientAndTorque[1][0].E(Double.NaN);
            gradientAndTorque[0][1].E(Double.NaN);
            gradientAndTorque[1][1].E(Double.NaN);
            return gradientAndTorque;
        }
        Vector site0 = space.makeVector();
        Vector site1 = space.makeVector();
        Vector drij = space.makeVector();
        Vector torque = space.makeVector();
        gradientAndTorque[0][0].E(0);
        gradientAndTorque[0][1].E(0);
        gradientAndTorque[1][0].E(0);
        gradientAndTorque[1][1].E(0);
        for (int i=0; i<5; i++) {
            int ii = siteID[i];
            site0.Ea1Tv1(pos[i], or0);
            site0.PE(cm0);
            for (int j=0; j<5; j++) {
                int jj = siteID[j];
                site1.Ea1Tv1(pos[j], or1);
                site1.PE(cm1);
                drij.Ev1Mv2(site1, site0);
                double rij2 = drij.squared();            
                double rij = Math.sqrt(rij2);                    
                double ar = alpha[ii][jj]*rij;                
                double rduExpdr = -A[ii][jj]*ar*Math.exp(-ar);

                double sum = 1;
                double dsum = 0;
                double br = b[ii][jj]*rij;
                double term = 1;
                double dterm = 1;                
                for (int k=1; k<=6; k++) {
                    term *= br/k;
                    if (k==1) {
                        dterm = br;
                    }
                    else {
                        dterm *= br/(k-1);
                    }                    
                    sum += term;
                    dsum += dterm;                    
                }
                if (sum==1) {
                    gradientAndTorque[0][0].E(Double.NaN);
                    gradientAndTorque[1][0].E(Double.NaN);
                    gradientAndTorque[0][1].E(Double.NaN);
                    gradientAndTorque[1][1].E(Double.NaN);
                    return gradientAndTorque;
                }                
                double rij6 = rij2*rij2*rij2;
                double expbr = Math.exp(-br);
                double rdu6dr = (-expbr*br*sum + expbr*dsum + 6*(1-expbr*sum))*c6[ii][jj]/rij6;                
                double rduChargedr = -q[ii]*q[jj]/rij;                
                double rdudr = Kelvin.UNIT.toSim(rduExpdr + rdu6dr + rduChargedr);
                
                // we're done with drij.  replace it with gradient on 1 (force on 0)
//                System.out.println("g drij "+drij+" "+rdudr/rij2);
                drij.TE(rdudr/rij2);
                gradientAndTorque[0][1].PE(drij);
//                System.out.println("grad 1 "+gradientAndTorque[0][1]);
                gradientAndTorque[0][0].ME(drij);
                torque.Ea1Tv1(pos[j], or1);
                torque.XE(drij);
                gradientAndTorque[1][1].ME(torque);
                torque.Ea1Tv1(pos[i], or0);
                torque.XE(drij);
                gradientAndTorque[1][0].PE(torque);
                if (drij.isNaN() || gradientAndTorque[0][1].isNaN()) {
                    throw new RuntimeException();
                }
            }
        }
        return gradientAndTorque;
    }


    public double virial(IAtomList atoms) {
        return 0;
    }


    public Vector[] gradient(IAtomList atoms) {
        return gradientAndTorque(atoms)[0];
    }


    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradientAndTorque(atoms)[0];
    }
    
}

