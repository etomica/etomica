package simulate;
import simulate.units.*;

/**
 * This potential is generally designed to work with molecules of three atoms.  The first
 * molecule acts as an attractive site, the second has a square-well
 * potential with other molecules, and the rest have ideal gas interactions.  The third atom
 * interacts with the others usually with a tether potential, and it serves to keep the attractive site from
 * moving through the center of the molecule.  
 *
 * @author Rob Riggleman
 */

public class P2Association extends Potential2 {
    
    private final PotentialIdealGas potentialIdealGas;
    private final PotentialSquareWell potentialSquareWell;
    private final PotentialAssociation potentialAssociation;
    private double lambda = 1.5, epsilon, coreDiameter, siteWellDiameter = 4.5;
    
    public P2Association() {
        this(Simulation.instance);
    }
    public P2Association(Simulation sim) {
        super(sim);
        setCoreDiameter(Default.ATOM_SIZE);
        setEpsilon(Default.POTENTIAL_WELL);
        setLambda(lambda);
        setSiteWellDiameter(siteWellDiameter);
        potentialIdealGas = new PotentialIdealGas();
        potentialSquareWell = new PotentialSquareWell(coreDiameter, lambda, epsilon);
        potentialAssociation = new PotentialAssociation(siteWellDiameter, epsilon);
    }
    
    public final Potential getPotential(Atom a1, Atom a2) {
        if (a1.atomIndex() == 1 && a2.atomIndex() == 1) {
            return potentialSquareWell;
        }
        else if (a1.atomIndex() == 0 && a2.atomIndex() == 0) {
            return potentialAssociation;
        }
        else {return potentialIdealGas;}
    }
    
    /**
     * Accessor method for the core diameter of the square-well potential for the central atoms.
     */
    public final void setCoreDiameter(double d) {
        coreDiameter = d;
        setPotentialCutoff(coreDiameter*lambda);
        potentialSquareWell.setCoreDiameter(coreDiameter);
    }
    /**
     * Accessor method for the core diameter of the square-well potential for the central atoms.
     */
    public final double getCoreDiameter() {return coreDiameter;}
    public Dimension getCoreDiameterDimension() {return Dimension.LENGTH;}
    
    /**
     * Accessor method for the well-width multiplier of the square-well potential for the central atoms.
     */
    public double getLambda() {return lambda;}
    /**
     * Accessor method for the well-width multiplier of the square-well potential for the central atoms.
     */
    public void setLambda(double lam) {
        lambda = lam;
        if (coreDiameter*lambda > siteWellDiameter) setPotentialCutoff(coreDiameter*lambda);
        if (potentialSquareWell != null) potentialSquareWell.setLambda(lambda);
    }
    /**
     * Indicates that the well-width multiplier is dimensionless
     */
    public Dimension getLambdaDimension() {return Dimension.NULL;}
    
    /**
     * Accessor method for the well-width diameter of the attractive site on the molecule.
     */
    public double getSiteWellDiameter() {return siteWellDiameter;}
    /**
     * Accessor method for the well-width diameter of the attractive site on the molecule.
     */
    public void setSiteWellDiameter(double well) {
        siteWellDiameter = well;
        if (siteWellDiameter > coreDiameter*lambda) setPotentialCutoff(siteWellDiameter);
        potentialAssociation.setWellDiameter(siteWellDiameter);
    }
    /**
     * Indicates that the attractive well diameter has dimension of length
     */
    public Dimension getWellDiameterDimension() {return Dimension.LENGTH;}
    
    /**
     * Accessor method for the well depth of both the central square well and the attractive site
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Indicates that the well depth has dimensions of energy
     */
    public Dimension getEpsilonDimension() {return Dimension.ENERGY;}
    /**
     * Accessor method for the well depth of both the central square well and the attractive site
     */
    public void setEpsilon(double e) {
        epsilon = e;
        potentialSquareWell.setEpsilon(epsilon);
        potentialAssociation.setEpsilon(epsilon);
    }
}