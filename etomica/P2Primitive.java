package etomica;

import etomica.electrostatics.*;

/**
 * Potential for the primitive model of an electrolyte.
 * The primitive model is a simple Coulombic charge with a hard-sphere repulsion.
 * Charges are ascribed to the Atoms via their AtomType.electroType field
 * Hard-sphere diameter is defined by this potential
 */
public class P2Primitive extends Potential2 implements EtomicaElement {

    public String getVersion() {return "P2Primitive:01.07.08/"+Potential2SoftSpherical.VERSION;}

    private double sigma, sigmaSquared;
    private double cutoffRadius, cutoffRadiusSquared;
    private double cutoff;
    private double eLRC, pLRC;  //multipliers for long-range correction to energy and pressure, resp.
    private Space.Vector force;

    public P2Primitive() {
        this(Simulation.instance.hamiltonian.potential, Default.ATOM_SIZE);
    }
    public P2Primitive(double sigma) {
        this(Simulation.instance.hamiltonian.potential, sigma);
    }
    public P2Primitive(PotentialGroup parent, double sigma) {
        super(parent);
        setSigma(sigma);
        force = parentSimulation().space().makeVector();
    }
 
   /**
    * Returns primitive-model energy.
    * Return infinity if overlap is true, and zero if separation is greater than cutoff.
    */
    public double energy(AtomPair pair) {
        double r2 = pair.r2();
        if(r2 < sigmaSquared) {return Double.MAX_VALUE;}
        else {
            double z1 = ((Monopole)pair.atom1.type.electroType()).getZ();
            double z2 = ((Monopole)pair.atom2.type.electroType()).getZ();
//            return z1*z2/Math.sqrt(r2);
            return -0.5*z1*z2*Math.log(r2);
        }
    }
    /** 
     * Force that atom2 exerts on atom1
     * This method has not been checked and may be incorrect
     */
    public Space.Vector force(AtomPair pair) {
        double r2 = pair.r2();
        if(r2 > cutoffRadiusSquared) {force.E(0.0);}
        else if(r2 < sigmaSquared) {force.E(Double.MAX_VALUE);}
        else {
            double z1 = ((Monopole)pair.atom1.type.electroType()).getZ();
            double z2 = ((Monopole)pair.atom2.type.electroType()).getZ();
            double c = z1*z2/Math.sqrt(r2);
            force.E(pair.dr());
            force.TE(-c/r2);
        }
        return force;
    }   
    
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        if( !(pc instanceof Potential2Calculation) ) return;
        iterator.reset(id);
        ((Potential2Calculation)pc).calculate(iterator, this); 
    }//end of calculate
     
    /**
     * Accessor method for the size of the repulsive core of the primitive model
     */
    public double getSigma() {return sigma;}
    /**
     * Accessor method for the size of the repulsive core of the primitive model
     */
    public final void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
    }
    public etomica.units.Dimension getSigmaDimension() {return etomica.units.Dimension.LENGTH;}

}
  