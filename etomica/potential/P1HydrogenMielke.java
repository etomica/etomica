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
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.units.Kelvin;

public class P1HydrogenMielke implements IPotential{
//    1 Eh (hartree) = 27.2113961eV = 627.5096 kcal/mol = 219474.7 cm^-1 

    protected IBoundary boundary;
    protected final double c6 = -6.499027;
    protected final double c8 = -124.3991;
    protected final double c10 = -3285.828;
    protected final double beta = 0.2;
    protected final double re = 1.401;
    protected final double alpha = 3.980917850296971;
    protected final double[] xlp = {-1.709663318869429E-1,-6.675286769482015E-1,-1.101876055072129E+0,-1.106460095658282E+0,-7.414724284525231E-1,-3.487882731923523E-1,-1.276736255756598E-1,-5.875965709151867E-2,-4.030128017933840E-2,-2.038653237221007E-2,-7.639198558462706E-4,2.912954920483885E-3,-2.628273116815280E-4,-4.622088855684211E-4,1.278468948126147E-4,-1.157434070240206E-5,-2.609691840882097E-12};
    
    public static void main(String[] args) {
        ISpace space = Space3D.getInstance();        
        P1HydrogenMielke pot1 = new P1HydrogenMielke(space);
        int maxi = 1000;
//        System.out.println("1hartree to sim = "+Hartree.UNIT.toSim(1)+", 1K to sim = "+Kelvin.UNIT.toSim(1));

        try {
            FileWriter file = new FileWriter("MielkeFit1.dat");
            
            for (int i=0; i<maxi; i++) {
                double r1 = 0.2+1.8*i/maxi;
                double E1 = Hartree.UNIT.fromSim(pot1.vH2(r1))*627.5096;
                double E2 = Kelvin.UNIT.fromSim(pot1.vH2(r1));
                file.write(r1+" "+E1+" "+E2+"\n");
            }
            file.close();
            
        } catch(IOException e) {
            throw new RuntimeException("caught IOException: " + e.getMessage());            
        }       
        
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

    public double vH2(double pr) {
        if (false) return 100000*(pr-1)*(pr-1);
        // H2 singlet potential for FCI/CBS
        double r = BohrRadius.UNIT.fromSim(pr);
        double damp=1.0-Math.exp(-beta*r*r);
        double xd = damp/r ;
        double xd2 = xd*xd;
        double xd6 = xd2*xd2*xd2;
        double xd8 = xd6*xd2;
        double xd10 = xd8*xd2;
        double prefac = Math.exp(- alpha*(r-re));
        double eLong = c6*xd6 + c8*xd8 + c10*xd10;
        double eShort = 0.00;        
        for (int i=0; i<17; i++) {
            eShort += xlp[i]*prefac;
            prefac *= (r-re);
        }
        double f = Hartree.UNIT.toSim(eLong + eShort);
        return f;
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
    
}
