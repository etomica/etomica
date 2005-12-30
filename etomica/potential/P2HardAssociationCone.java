package etomica.potential;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.CoordinatePair;
import etomica.space.ICoordinateAngular;
import etomica.space.Space;
import etomica.space.Vector;
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
    private double wellcutoffFactor;
    private double /*wellCutoff, */wellCutoffSquared;
    private double sigma, sigmaSquared;
    private double epsilon, epsilon4, wellEpsilon;
    private double cutoffLJSquared, cutoffFactor;
    private Vector e1;
    private Vector e2;
    private double theta, ec2;
    private final CoordinatePair cPair;
    
    public P2HardAssociationCone(Simulation sim) {
        this(sim.space, sim.getDefaults().atomSize, sim.getDefaults().potentialWell, 
                sim.getDefaults().potentialCutoffFactor);
    }
    
    public P2HardAssociationCone(Space space, double sigma, double epsilon, double cutoffFactorLJ) {
        super(space);
        e1 = space.makeVector();
        e2 = space.makeVector();
        cPair = new CoordinatePair(space);

        setSigma(sigma);
        setEpsilon(epsilon);
        setCutoffFactorLJ(cutoffFactorLJ);
//        setWellCutoff(getSigma());
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
    public double energy(AtomSet pair) {
    	    cPair.reset((AtomPair)pair);
        double eTot = 0.0;
        double r2 = cPair.r2();
                 
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
            ((ICoordinateAngular)((AtomPair)pair).atom0.coord).orientation().convertToSpaceFrame(e1);
            double er1 = e1.dot(cPair.dr());
                       
            if ( er1 > 0.0 && er1*er1 > ec2*r2) {
                e2.E(0.);
                e2.setX(0,1);
                ((ICoordinateAngular)((AtomPair)pair).atom1.coord).orientation().convertToSpaceFrame(e2);
                double er2 = e2.dot(cPair.dr());
                if(er2 < 0.0 && er2*er2 > ec2*r2) eTot -= wellEpsilon;
            }
        }
        return eTot;
    }//end of energy
    
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
       // calculateLRC();
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
        wellcutoffFactor=wcut;
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
    * Accessor method for attractive-well diameter.
    */
//    public double getWellCutoff() {return wellCutoff;}
    /**
    * Accessor method for attractive-well diameter.
    */
/*    public void setWellCutoff(double wcut) {
        wellCutoff = wcut;
        wellCutoffSquared = wcut*wcut;
    }
 */         
    public static final Dimension getWellCutoffDimension() {return Length.DIMENSION;}
    
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
        ec2    = Math.cos(theta);
        ec2   = ec2*ec2;
    }
    public Dimension getThetaDimension() {return Angle.DIMENSION;}

    public void setPhase(Phase phase) {
        cPair.setNearestImageTransformer(phase.getBoundary());
    }

    /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
	    IntegratorMC integrator1 = new IntegratorMC();
	    SpeciesSpheresRotating speciesSpheres1 = new SpeciesSpheresRotating(20);
	    Phase phase1 = new Phase();
	    PotentialAssociationCone potential = new PotentialAssociationCone();
	//    potential.setEpsilon(0.0);
	    potential.setWellCutoff(1.5*potential.getSigma());
	    P2SimpleWrapper p2StrongFluids1 = new P2SimpleWrapper(potential);
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    MeterEnergy meterEnergy = new MeterEnergy();
	    meterEnergy.setPhase(phase1);
	    DisplayBox displayEnergy = new DisplayBox();
	    displayEnergy.setMeter(meterEnergy);
	    displayEnergy.setUnit(new etomica.units.Unit(etomica.units.LennardJones.Energy.UNIT));
	    MCMove mcAtom = new MCMoveAtom();
	    MCMove mcRotate = new MCMoveRotate();
	    integrator1.add(mcAtom);
	    integrator1.add(mcRotate);
		Simulation.instance.setBackground(java.awt.Color.yellow);
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
    
        */
}