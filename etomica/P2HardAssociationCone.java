package etomica;
import etomica.units.Dimension;

/**
 * Lennard-Jones potential with a square-well cone of attraction. 
 *
 * @author Jayant K. Singh
 */

public class P2HardAssociationCone extends Potential2 implements EtomicaElement {
    private double wellCutoff, wellCutoffSquared;
    private double sigma, sigmaSquared;
    private double epsilon, epsilon4, wellEpsilon;
    private double cutoffLJSquared, cutoffFactor;
    private Space.Vector e1;
    private Space.Vector e2;
    private double theta, ec2;
    
    public P2HardAssociationCone() {this(Simulation.instance.hamiltonian.potential);}
    
    public P2HardAssociationCone(PotentialGroup parent) {
        super(parent);
        e1 = parentSimulation().space().makeVector();
        e2 = parentSimulation().space().makeVector();

        setSigma(Default.ATOM_SIZE);
        setEpsilon(Default.POTENTIAL_WELL);
        setCutoffFactorLJ(Default.POTENTIAL_CUTOFF_FACTOR);
        setWellCutoff(getSigma());
        setWellEpsilon(8.0*getEpsilon());
        setTheta(etomica.units.Degree.UNIT.toSim(27.0));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Lennard-Jones core with an anisotropic, cone-shaped region of square-well attraction");
        return info;
    }


    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        if( !(pc instanceof Potential2Calculation) ) return;
        iterator.reset(id);
        ((Potential2Calculation)pc).calculate(iterator, this); 
    }//end of calculate

 /**
  * Returns the pair potential energy.
  */
    public double energy(AtomPair pair) {
        double eTot = 0.0;
        Atom a1 = pair.atom1();
        Atom a2 = pair.atom2();
        double r2 = pair.r2();
                 
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
            e1.setComponent(0,1);
            ((Space.Coordinate.Angular)a1.coord).orientation().convertToSpaceFrame(e1);
            double er1 = e1.dot(pair.dr());
                       
            if ( er1 > 0.0 && er1*er1 > ec2*r2) {
                e2.E(0.);
                e2.setComponent(0,1);
                ((Space.Coordinate.Angular)a2.coord).orientation().convertToSpaceFrame(e2);
                double er2 = e2.dot(pair.dr());
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
    public static final Dimension getSigmaDimension() {return Dimension.LENGTH;}

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
    public static final Dimension getCutoffFactorLJDimension() {return Dimension.NULL;}
   
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
    public static final Dimension getEpsilonDimension() {return Dimension.ENERGY;}
    
    /**
    * Accessor method for attractive-well diameter.
    */
    public double getWellCutoff() {return wellCutoff;}
    /**
    * Accessor method for attractive-well diameter.
    */
    public void setWellCutoff(double wcut) {
        wellCutoff = wcut;
        wellCutoffSquared = wcut*wcut;
    }
          
    public static final Dimension getWellCutoffDimension() {return Dimension.LENGTH;}
    
    /**
    * Accessor method for attractive-well depth parameter.
    */
    public double getWellEpsilon() {return wellEpsilon;}
    /**
    * Accessor method for attractive-well depth parameter.
    */
    public void setWellEpsilon(double weps) {wellEpsilon = weps;}
          
    public static final Dimension getWellEpsilonDimension() {return Dimension.ENERGY;}
    
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
    public Dimension getThetaDimension() {return Dimension.ANGLE;}

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