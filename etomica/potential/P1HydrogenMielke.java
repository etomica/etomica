package etomica.potential;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.api.IVector;
import etomica.atom.AtomHydrogen;
import etomica.atom.IAtomOriented;
import etomica.chem.elements.Hydrogen;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.units.Kelvin;

public class P1HydrogenMielke implements IPotential{
//    1 Eh (hartree) = 27.2113961eV = 627.5096 kcal/mol = 219474.7 cm^-1 

    protected IBoundary boundary;
    protected final static double c6 = -6.499027;
    protected final static double c8 = -124.3991;
    protected final static double c10 = -3285.828;
    protected final static double beta = 0.2;
    protected final static double re = 1.401;
    protected final static double alpha = 3.980917850296971;
    protected final static double[] xlp = {-1.709663318869429E-1,-6.675286769482015E-1,-1.101876055072129E+0,-1.106460095658282E+0,-7.414724284525231E-1,-3.487882731923523E-1,-1.276736255756598E-1,-5.875965709151867E-2,-4.030128017933840E-2,-2.038653237221007E-2,-7.639198558462706E-4,2.912954920483885E-3,-2.628273116815280E-4,-4.622088855684211E-4,1.278468948126147E-4,-1.157434070240206E-5,-2.609691840882097E-12};
    
    public static void main(String[] args) {        
        int maxi = 1000;
//        System.out.println("1hartree to sim = "+Hartree.UNIT.toSim(1)+", 1K to sim = "+Kelvin.UNIT.toSim(1));

//        try {
//            FileWriter file = new FileWriter("MielkeFit1.dat");
//            
//            for (int i=0; i<maxi; i++) {
//                double r1 = 0.2+1.8*i/maxi;
//                double E1 = Hartree.UNIT.fromSim(pot1.vH2(r1))*627.5096;
//                double E2 = Kelvin.UNIT.fromSim(pot1.vH2(r1));
//                file.write(r1+" "+E1+" "+E2+"\n");
//            }
//            file.close();
//            
//        } catch(IOException e) {
//            throw new RuntimeException("caught IOException: " + e.getMessage());            
//        }
        double P = 256.0;
        double t = 12.471707671342523;
        double r1 = 0.741419798250176*1.001;
        double r0 = 0.741419798250176;
        double b1 = P1HydrogenMielke.vH2(r1);
        double b0 = P1HydrogenMielke.vH2(r0);
        System.out.println(b1+" "+b0);
        double lhs = b1 - b0;
        double rhs = (r1 - r0)*P1HydrogenMielke.dvdr(r0) + 0.5*(r1 - r0)*(r1 - r0)*P1HydrogenMielke.d2vdr2(r0);
        System.out.println(lhs+" "+rhs+" "+P1HydrogenMielke.dvdr(r0)+" "+P1HydrogenMielke.d2vdr2(r0));
        
//        for (int i=0; i<maxi; i++) {
//            double r = 0.2+9.8*i/maxi;
//            double E1 = Hartree.UNIT.fromSim(pot1.vH2(r1))*627.5096;
//            double E2 = Kelvin.UNIT.fromSim(pot1.vH2(r1));
//            System.out.println(r+" "+P1HydrogenMielke.vH2(r)+" "+P1HydrogenMielke.dvdrNew(r)+" "+P1HydrogenMielke.d2vdr2New(r));
//            System.out.println(r1+" "+pot1.vH2(r1));
//            file.write(r1+" "+E1+" "+E2+"\n");
//        }
        
    }

    public P1HydrogenMielke(ISpace space) {        
    }

    public double getRange() {    
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(IBox box) {    
        boundary = box.getBoundary();
    }

    public int nBody() {    
        return 1;
    }    

    public static double vH2(double pr) {       
        // H2 singlet potential for FCI/CBS
    	if (pr < 0) return Double.POSITIVE_INFINITY;
        double r = BohrRadius.UNIT.fromSim(pr);
        double damp=1.0-Math.exp(-beta*r*r);
        double xd = damp/r ;
        double xd2 = xd*xd;
        double xd6 = xd2*xd2*xd2;
        double xd8 = xd6*xd2;
        double xd10 = xd8*xd2;
        double prefac = Math.exp(- alpha*(r-re));
        double vLong = c6*xd6 + c8*xd8 + c10*xd10;
        double vShort = 0.00;        
        for (int i=0; i<17; i++) {
            vShort += xlp[i]*prefac;
            prefac *= (r-re);
        }
        double v = Hartree.UNIT.toSim(vLong + vShort);
        return v;
    }
    
    public static double dvdr (double pr) {
    	double r = BohrRadius.UNIT.fromSim(pr);
    	double r2 = r*r;
    	double z = beta*r2;
        double y = Math.exp(-z);
        double inverseRm = (1.0 - y)/r ;
        double inverseRm2 = inverseRm*inverseRm;
        double inverseRm5 = inverseRm2*inverseRm2*inverseRm;
        double inverseRm7 = inverseRm5*inverseRm2;
        double inverseRm9 = inverseRm7*inverseRm2;        
        double dinverseRmdr = (2*z*y + y - 1)/r2;
        double dvLongdr = dinverseRmdr*(6*c6*inverseRm5 + 8*c8*inverseRm7 + 10*c10*inverseRm9);
        double dvShortdr = -alpha*xlp[0];
        double fac = 1.0;
        for (int i=1; i<17; i++) {
            dvShortdr += xlp[i]*fac*(i - alpha*(r-re));
            fac *= (r-re);
        }
        dvShortdr *= Math.exp(-alpha*(r-re));
        double dv = Hartree.UNIT.toSim(dvLongdr + dvShortdr)/BohrRadius.UNIT.toSim(1);         
        return dv;
    }   
    
    public static double d2vdr2(double pr) {        
        double r = BohrRadius.UNIT.fromSim(pr);
        double r2 = r*r;
    	double z = beta*r2;
        double y = Math.exp(-z);        
        double inverseRm = (1.0 - y)/r ;
        double inverseRm2 = inverseRm*inverseRm;
        double inverseRm4 = inverseRm2*inverseRm2;
        double inverseRm5 = inverseRm4*inverseRm;
        double inverseRm6 = inverseRm5*inverseRm;
        double inverseRm7 = inverseRm6*inverseRm;
        double inverseRm8 = inverseRm7*inverseRm;
        double inverseRm9 = inverseRm8*inverseRm;
        double dinverseRmdr = (2*z*y + y - 1)/r2;
        double d2inverseRmdr2 = (2 - 2*y - 2*z*y - 4*z*z*y)/(r2*r);        
        double d2vL = d2inverseRmdr2*(6*c6*inverseRm5 + 8*c8*inverseRm7 + 10*c10*inverseRm9) + dinverseRmdr*dinverseRmdr*(30*c6*inverseRm4 + 56*c8*inverseRm6 + 90*c10*inverseRm8);
        double d2vS = alpha*alpha*xlp[0] + (alpha*alpha*(r-re) - 2*alpha)*xlp[1];
        double fac = 1;
        for (int i=2; i<17; i++) {
            d2vS += xlp[i]*fac*(alpha*alpha*(r-re)*(r-re) + i*(i-1) - 2*i*alpha*(r-re));
            fac *= (r-re);
        }
        d2vS *= Math.exp(-alpha*(r-re));
        double d2v = Hartree.UNIT.toSim(d2vL + d2vS)/(BohrRadius.UNIT.toSim(1)*BohrRadius.UNIT.toSim(1));
        return d2v;
    }
    
    
    public static class P1HydrogenMielkeAtomic extends P1HydrogenMielke implements IPotentialAtomic {
        public P1HydrogenMielkeAtomic(ISpace space) {
            super(space);     
        }

        public double energy(IAtomList atoms) {
            AtomHydrogen m0 = (AtomHydrogen)atoms.getAtom(0);
            double bL = m0.getBondLength();
            double f = vH2(bL);
            return f;
        } 
    }
    
    public static class P2HydrogenMielkeAtomic extends P1HydrogenMielke implements IPotentialAtomic {
    	public P2HydrogenMielkeAtomic(ISpace space) {
    		super(space);
    	}
    	
    	public int nBody(){
    		return 2;
    	}
    	
    	public double energy(IAtomList atoms) {
    		IAtom a0 = atoms.getAtom(0);
    		IAtom a1 = atoms.getAtom(1);
            double bL = Math.sqrt(a0.getPosition().Mv1Squared(a1.getPosition()));
            double f = vH2(bL);
            return f;
    	}
    }
    
}
