package simulate;
import java.io.*;
import java.beans.*;
import java.awt.*;

public class PotentialSoftDisk extends Potential {
    
    private int n, halfN;
    private double epsilon, sigma;
    private double neighborCutoff, neighborCutoffSquared;
    private double potentialCutoff, potentialCutoffSquared;
    private double skinThickness;
    
    
    public PotentialSoftDisk() {
        super();
        n = 12;
        halfN = n/2;
        setEpsilon(1.0);
        setSigma(0.1);
        setPotentialCutoff(0.9*Double.MAX_VALUE);  //no cutoff for potential
        setSkinThickness(3.0*sigma);
    }
    
    public boolean isNeighbor(SpeciesElement e1, SpeciesElement e2){
        MoleculeAtomic disk1 = (MoleculeAtomic)e1;
        MoleculeAtomic disk2 = (MoleculeAtomic)e2;        
        space.uEr1Mr2(r12,disk2.r,disk1.r);  //use instance method   //r2-r1
        return (Space.v1Dv2(r12,r12) < neighborCutoffSquared);  
    }
    
/*    
    public int getN() {return n;} // n = 12 for now
    public void setN(int n) {this.n = n;}
*/    
    public double[] force(SpeciesElement element1, SpeciesElement element2) {
      MoleculeAtomic disk1 = (MoleculeAtomic)element1;
      MoleculeAtomic disk2 = (MoleculeAtomic)element2;        
      space.uEr1Mr2(r12,disk2.r,disk1.r);  //use instance method   //r2-r1
      double r2 = Space.v1Dv2(r12,r12);
      double factor = sigma*sigma/r2;
      factor = factor*factor;  // (sig/r)^4
      factor = factor*factor*factor; //(sig/r)^12
      factor = factor*n*epsilon/r2;
      Space.uEa1Tv1(f, factor, r12);
      return f;
  }
  
  public double dfdr(double r2) {
    double factor = sigma*sigma/r2;
    factor = factor*factor;  // (sig/r)^4
    factor = factor*factor*factor; //(sig/r)^12
    return factor*n*(n+1)*epsilon/r2;
  }

  public PairInteraction computePairInteraction(SpeciesElement element1, SpeciesElement element2) {
      MoleculeAtomic disk1 = (MoleculeAtomic)element1;
      MoleculeAtomic disk2 = (MoleculeAtomic)element2;        
      space.uEr1Mr2(pair.rij,disk2.r,disk1.r);
      pair.rSquared = Space.v1Dv2(pair.rij, pair.rij);
      if(pair.rSquared < potentialCutoffSquared) {
        double factor = sigma*sigma/pair.rSquared;
        factor = factor*factor;  // (sig/r)^4
        factor = factor*factor*factor; //(sig/r)^12
        pair.energy = epsilon*factor;
        pair.virial = -n*pair.energy;
        factor = factor*n*epsilon/pair.rSquared;
        Space.uEa1Tv1(pair.force, factor, pair.rij);
      }
      else {
        pair.energy = pair.virial = 0.0;
        Space.uEa1(pair.force, 0.0);
      }
      return pair;
  }
 /*   
  public void setNeighborShell(double nbr) {
    neighborShell = nbr;
    neighborCutoffSquared = neighborShell*neighborShell*sigma*sigma;
  }
  public double getNeighborShell() {return neighborShell;}
  */
  
  public void setPotentialCutoff(double cutoff) {
    potentialCutoff = cutoff;
    potentialCutoffSquared = cutoff*cutoff;
    neighborCutoff = potentialCutoff + skinThickness;
    neighborCutoffSquared = neighborCutoff*neighborCutoff;
  }
  public double getPotentialCutoff() {return potentialCutoff;}
  
  public void setSkinThickness(double thickness) {
    skinThickness = thickness;
    neighborCutoff = potentialCutoff + skinThickness;
    neighborCutoffSquared = neighborCutoff*neighborCutoff;
  }
  public double getSkinThickness() {return skinThickness;}
  
  public void setEpsilon(double eps) {epsilon = eps;}
  public double getEpsilon() {return epsilon;}
  
  public void setSigma(double sig) {sigma = sig;}
  public double getSigma() {return sigma;}
    
    
}