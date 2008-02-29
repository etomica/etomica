package etomica.potential;

import etomica.EtomicaInfo;
import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;

/**
 * Hard core with an attractive tail that goes to zero linearly with r.
 * 
 * @author Jhumpa Adhikari
 */

public class P2TriangleWell extends Potential2 {

    public P2TriangleWell(Space space) {
        this(space, 1.0, 2.0, 1.0);
    }
  
    public P2TriangleWell(Space space, double coreDiameter, double lambda, double epsilon) {
        super(space);
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilon(epsilon);
        force = space.makeVector();
        dr = space.makeVector();
    }
  
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard core with a surrounding region of constant-force attraction");
        return info;
    }

    public double energy(AtomSet pair) {
        IAtomPositioned atom0 = (IAtomPositioned)pair.getAtom(0);
        IAtomPositioned atom1 = (IAtomPositioned)pair.getAtom(1);
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        nearestImageTransformer.nearestImage(dr);

        double r2 = dr.squared();
       
        if(r2 < coreDiameterSquared)
            return Double.POSITIVE_INFINITY;

        if(r2 > wellDiameterSquared)
            return 0.0;
            
        double r1 = Math.sqrt(r2);
            
        return (epsilon/(lambda - 1.0))*((r1/coreDiameter)- lambda);
    }
 

    // what could call this?
    public IVector force(AtomSet pair){
        
        IAtomPositioned atom0 = (IAtomPositioned)pair.getAtom(0);
        IAtomPositioned atom1 = (IAtomPositioned)pair.getAtom(1);
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        nearestImageTransformer.nearestImage(dr);

        double r2 = dr.squared();
        if(r2 > wellDiameterSquared){
            force.E(0.0);
        }
        else {
            force.E(dr);
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

    public void setBox(Box box) {
        nearestImageTransformer = box.getBoundary();
    }

    private static final long serialVersionUID = 1L;
    private double coreDiameter, coreDiameterSquared;
    private double wellDiameter, wellDiameterSquared;
    private double lambda; //wellDiameter = coreDiameter * lambda ;lambda is well width
    private double epsilon;
    private double constant;
    private final IVector force;
    private final IVector dr;
    private NearestImageTransformer nearestImageTransformer;
}

  