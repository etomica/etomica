package simulate;

public class PotentialLJ extends Potential implements PotentialSoft {

  private double sigma, sigmaSquared;
  private double cutoffDiameter, cutoffDiameterSquared;
  private double epsilon, cutoff;
  private double epsilon4, epsilon24;

  public PotentialLJ(double sigma, double epsilon, double cutoff) {
    setSigma(sigma);
    setCutoff(cutoff);
    setEpsilon(epsilon);
  }
   
    public double energy(PhaseSpace.AtomPair pair) {
        double r2 = pair.r2();
        if(r2 > cutoffDiameterSquared) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return epsilon4*s6*(s6 - 1.0);
        }
    }
 
    public double virial(PhaseSpace.AtomPair pair) {  //not carefully checked for correctness
        double r2 = pair.r2();
        if(r2 > cutoffDiameterSquared) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return -epsilon24*s6*(2.0*s6 - 1.0);
        }
    }
    
/*    public double[] force(PhaseSpace.AtomPair pair) {  //not carefully checked for correctness
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
    }*/
        
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
  