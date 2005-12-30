package etomica.potential;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.CoordinatePair;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Hard core with an attractive tail that goes to zero linearly with r.
 * 
 * @author Jhumpa Adhikari
 */

public class P2TriangleWell extends Potential2 implements EtomicaElement {

  private double coreDiameter, coreDiameterSquared;
  private double wellDiameter, wellDiameterSquared;
  private double lambda; //wellDiameter = coreDiameter * lambda ;lambda is well width
  private double epsilon;
  private double constant;
  private final Vector force;
protected CoordinatePair cPair;

  public P2TriangleWell(Simulation sim) {
    this(sim.space, sim.getDefaults().atomSize, 
            sim.getDefaults().potentialCutoffFactor, sim.getDefaults().potentialWell);
  }
  
  public P2TriangleWell(Space space, double coreDiameter, double lambda, double epsilon) {
    super(space);
    setCoreDiameter(coreDiameter);
    setLambda(lambda);
    setEpsilon(epsilon);
    force = space.makeVector();
  }
  
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard core with a surrounding region of constant-force attraction");
        return info;
    }

    public double energy(AtomSet pair) {

    	cPair.reset((AtomPair)pair);
        double r2 = cPair.r2();
       
        if(r2 < coreDiameterSquared)
            return Double.POSITIVE_INFINITY;

        if(r2 > wellDiameterSquared)
            return 0.0;
            
        double r1 = Math.sqrt(r2);
            
        return (epsilon/(lambda - 1.0))*((r1/coreDiameter)- lambda);
    }
 

    public Vector force(AtomSet pair){
        
    	cPair.reset((AtomPair)pair);
        double r2 = cPair.r2();
        if(r2 > wellDiameterSquared){
            force.E(0.0);
        }
        if(r2 < wellDiameterSquared){
            force.E(cPair.dr());
            force.TE(constant/Math.sqrt(r2));//lambda > 1.0
            
        }
        return force;
    }
     
    public double getCoreDiameter() {return coreDiameter;}
    public final void setCoreDiameter(double c) {
        coreDiameter = c;
        coreDiameterSquared = c*c;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
        constant = epsilon/(coreDiameter*(1.0 - lambda));
    }
    public final etomica.units.Dimension getCoreDiameterDimension() {
        return etomica.units.Length.DIMENSION;
    }
    
    public double getRange() {
        return wellDiameter;
    }

    public double getLambda() {return lambda;}
    public final void setLambda(double lam) {
        lambda = lam;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
        constant = epsilon/(coreDiameter*(1.0 - lambda));
    }
    public final etomica.units.Dimension getLambdaDimension() {
        return etomica.units.Null.DIMENSION;
    }

    public double getEpsilon() {return epsilon;}
    public final void setEpsilon(double eps) {
        epsilon = eps;
        constant = epsilon/(coreDiameter*(1.0 - lambda));
    }
    public final etomica.units.Dimension getEpsilonDimension() {
        return etomica.units.Energy.DIMENSION;
    }

    public void setPhase(Phase phase) {
        cPair.setNearestImageTransformer(phase.getBoundary());
    }
}

  