package etomica.potential;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateAngular;
import etomica.space.IVectorRandom;
import etomica.space.IVector;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.units.Angle;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;
import etomica.units.Null;

/**
 * Lennard-Jones potential with a square-well cone of attraction. 
 *
 * @author Jayant K. Singh
 */

public class P2HardAssociationCone extends Potential2 implements EtomicaElement {
    private static final long serialVersionUID = 1L;
    private double wellcutoffFactor;
    private double wellCutoffSquared;
    private double sigma, sigmaSquared;
    private double epsilon, epsilon4, wellEpsilon;
    private double cutoffLJSquared, cutoffFactor;
    private IVectorRandom e1;
    private IVectorRandom e2;
    private double theta, ec2;
    private final IVector[] eArr = new IVector[1];
    private final IVectorRandom dr;
    private NearestImageTransformer nearestImageTransformer;
    
    public P2HardAssociationCone(Simulation sim) {
        this(sim.getSpace(), sim.getDefaults().atomSize, sim.getDefaults().potentialWell, 
                sim.getDefaults().potentialCutoffFactor);
    }
    
    public P2HardAssociationCone(Space space, double sigma, double epsilon, double cutoffFactorLJ) {
        super(space);
        e1 = space.makeVector();
        e2 = space.makeVector();
        dr = space.makeVector();

        setSigma(sigma);
        setEpsilon(epsilon);
        setCutoffFactorLJ(cutoffFactorLJ);
        setWellCutoffFactor(1.0);
        setWellEpsilon(8.0*getEpsilon());
        setTheta(etomica.units.Degree.UNIT.toSim(27.0));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Lennard-Jones core with an anisotropic, cone-shaped region of square-well attraction");
        return info;
    }
    
    /**
     * Returns infinity.
     */
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }


    /**
     * Returns the pair potential energy.
     */
    public double energy(AtomSet atoms) {
        AtomPair pair = (AtomPair)atoms;
        ICoordinateAngular coord0 = (ICoordinateAngular)((AtomLeaf)pair.atom0).getCoord();
        ICoordinateAngular coord1 = (ICoordinateAngular)((AtomLeaf)pair.atom1).getCoord();
        dr.Ev1Mv2(coord1.getPosition(),coord0.getPosition());
        nearestImageTransformer.nearestImage(dr);
        double r2 = dr.squared();
        double eTot = 0.0;
                 
        if(r2 > cutoffLJSquared) {
            eTot = 0.0;
        }
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            eTot = epsilon4*s6*(s6 - 1.0);
        }
                  
        if (r2 < wellCutoffSquared) {
            e1.E(0.);
            e1.setX(0,1);
            eArr[0] = e1;
            coord0.getOrientation().convertToSpaceFrame(eArr);
            double er1 = e1.dot(dr);

            if ( er1 > 0.0 && er1*er1 > ec2*r2) {
                e2.E(0.);
                e2.setX(0,1);
                eArr[0] = e2;
                coord1.getOrientation().convertToSpaceFrame(eArr);
                double er2 = e2.dot(dr);
                if(er2 < 0.0 && er2*er2 > ec2*r2) eTot -= wellEpsilon;
            }
        }
        return eTot;
    }
    
    /**
     * Accessor method for Lennard-Jones size parameter
     */
    public double getSigma() {return sigma;}
    /**
     * Accessor method for Lennard-Jones size parameter
     */
    public void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
        setCutoffFactorLJ(cutoffFactor);
    }
    public static final Dimension getSigmaDimension() {return Length.DIMENSION;}

    /**
    * Accessor method for Lennard-Jones cutoff distance; divided by sigma
    * @return cutoff distance, divided by size parameter (sigma)
    */
    public double getCutoffFactorLJ() {return cutoffFactor;}
    /**
     * Accessor method for Lennard-Jones cutoff distance; divided by sigma
     * @param rc cutoff distance, divided by size parameter (sigma)
     */
    public void setCutoffFactorLJ(double rc) {  
        cutoffFactor = rc;
        double cutoffLJ = sigma*cutoffFactor;
        cutoffLJSquared = cutoffLJ*cutoffLJ;
    }
    public static final Dimension getCutoffFactorLJDimension() {return Null.DIMENSION;}
   
    /**
    * Accessor method for attractive-well diameter divided by sigma
    */
    public double getWellCutoffFactor() {return wellcutoffFactor;}
    /**
    * Accessor method for attractive-well diameter divided by sigma;
    */
    public void setWellCutoffFactor(double wcut) {
        wellcutoffFactor = wcut;
        double wellCutoff = sigma*wcut;
        wellCutoffSquared = wellCutoff*wellCutoff;
    }
          
    public static final Dimension getWellCutoffFactorDimension() {return Null.DIMENSION;}

    /**
    * Accessor method for Lennard-Jones energy parameter
    */ 
    public double getEpsilon() {return epsilon;}
    /**
    * Accessor method for depth of well
    */
    public void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = 4.0 * eps;
    }
    public static final Dimension getEpsilonDimension() {return Energy.DIMENSION;}
    
    /**
    * Accessor method for attractive-well depth parameter.
    */
    public double getWellEpsilon() {return wellEpsilon;}
    /**
    * Accessor method for attractive-well depth parameter.
    */
    public void setWellEpsilon(double weps) {wellEpsilon = weps;}
          
    public static final Dimension getWellEpsilonDimension() {return Energy.DIMENSION;}
    
    /**
     * Accessor method for angle describing width of cone.
     */
    public double getTheta() {return theta;}
    
    /**
     * Accessor method for angle (in radians) describing width of cone.
     */
    public void setTheta(double t) {
        theta = t;
        ec2   = Math.cos(theta);
        ec2   = ec2*ec2;
    }
    public Dimension getThetaDimension() {return Angle.DIMENSION;}

    public void setPhase(Phase phase) {
        nearestImageTransformer = phase.getBoundary();
    }
}
