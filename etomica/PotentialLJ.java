package simulate;

import java.awt.*;
import java.beans.Beans;

public class PotentialLJ extends Potential implements PotentialSoft {

  private double sigma, sigmaSquared;
  private double cutoffDiameter, cutoffDiameterSquared;
  private double epsilon, cutoff;
  private double epsilon4, epsilon24;
  private transient final double[] r12 = new double[Space.D];
  private transient final double[] f12 = new double[Space.D];
//  private transient final double[] v12 = new double[Space.D];
//  private transient PairInteraction pair = new PairInteraction();

  public PotentialLJ(double sigma, double epsilon, double cutoff) {
    setSigma(sigma);
    setCutoff(cutoff);
    setEpsilon(epsilon);
  }

 /* public PairInteraction computePairInteraction(Molecule molecule1, Molecule molecule2) {
    MoleculeAtomic disk1 = (MoleculeAtomic)molecule1;
    MoleculeAtomic disk2 = (MoleculeAtomic)molecule2;        
    space.uEr1Mr2(pair.rij,disk2.r,disk1.r);
    pair.rSquared = Space.v1Dv2(pair.rij, pair.rij);
    pair.energy = (pair.rSquared < wellDiameterSquared) ? -epsilon : 0.0;
    return pair;
  }
 */
   
    public double energy(AtomC atom1, AtomC atom2) {
        double r2 = parentPhase.space.r1Mr2_S(atom1.r, atom2.r);
        if(r2 > cutoffDiameterSquared) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return epsilon4*s6*(s6 - 1.0);
        }
    }
 
    public double virial(AtomC atom1, AtomC atom2) {  //not carefully checked for correctness
        double r2 = parentPhase.space.r1Mr2_S(atom1.r, atom2.r);
        if(r2 > cutoffDiameterSquared) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return -epsilon24*s6*(2.0*s6 - 1.0);
        }
    }
    
    public double[] force(AtomC atom1, AtomC atom2) {  //not carefully checked for correctness
        parentPhase.space.uEr1Mr2(r12,atom1.r, atom2.r);
        double r2 = Space.v1S(r12);
        if(r2 > cutoffDiameterSquared) {Space.uEa1(f12,0.0);}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            double factor = -epsilon24*s6*(2.0*s6 - 1.0)/r2;
            Space.uEa1Tv1(f12,factor,r12);
        }
        return f12;
    }
        
    public double getSigma() {return sigma;}
    public void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
        setCutoff(cutoff);
    }

    public double getCutoff() {return cutoff;}
    public void setCutoff(double c) {
        cutoff = c;
        cutoffDiameter = sigma*cutoff;
        cutoffDiameterSquared = cutoffDiameter*cutoffDiameter;
    }
    
    public double getEpsilon() {return epsilon;}
    public void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = eps*4.0;
        epsilon24 = eps*24.0;
    }
}
  