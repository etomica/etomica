package simulate;

import java.awt.*;
import java.beans.Beans;

public class PotentialLennardJones extends Potential {

  private double sigma, sigmaSquared;
  private double epsilon, fourEpsilon;
  private double rCut, rCutSquared;
  private transient final double[] r12 = new double[Space.D];
  private transient PairInteraction pair = new PairInteraction();

  public PotentialLennardJones(double sigma, double epsilon) {
    setSigma(sigma);
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
  
    public double energy(Atom atom1, Atom atom2) {
        double r2 = sigmaSquared/space.r1Mr2_S(atom1.r, atom2.r);
        double r6 = r2*r2*r2;
        return fourEpsilon*r6*(r6-1.0);
    }
 

    public double getSigma() {return sigma;}
    public void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
    }
    
    public double getEpsilon() {return epsilon;}
    public void setEpsilon(double eps) {
        epsilon = eps;
        fourEpsilon = 4.0*epsilon;
    }

}
  