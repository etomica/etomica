package simulate; 
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.*;

/**
 * Gravity potential field, in which the force is proportional to the mass.
 * Direction and magnitude of gravity field is adjustable.  Default is value for earth
 * acting in the (D-1)-coordinate direction (e.g., the y-direction in 2D)
 */
public class PotentialFieldGravity extends PotentialField implements PotentialField.Soft {
  
    private final Space.Vector gVector;
    private final Space.Vector gForce;
    
    /**
     * Constructs gravitational potential field in "down" direction with magnitude equal to earth value
     */
    public PotentialFieldGravity(Phase p) {
        this(p, 1.0);
    }
    /**
     * Constructs gravity field in "down" direction, with magnitude given by multiplier of earth value.
     * Multiplier must be very large (of order 10^12) for field to have any noticable effect
     *
     * @param p phase in which this field is acting
     * @param multiplier magnitude of field as a multiple of value on Earth
     */
    public PotentialFieldGravity(Phase p, double multiplier) {
        super(p);
        gForce = p.parentSimulation().space().makeVector();
        gVector = p.parentSimulation().space().makeVector();
        gVector.setComponent(p.parentSimulation().space().D()-1, multiplier*Constants.G); //set last component of g-vector to g of earth
    }
    
    /**
     * Copies the elements of the given vector into the gravity vector
     */
    public void setGVector(Space.Vector g) {gVector.E(g);}
    /**
     * Accessor method for the gravity vector
     */
    public Space.Vector getGVector() {return gVector;}
    
    /**
     * Returns the energy due to the interaction of the atom with the field.
     */
    public double energy(Atom atom) {
        return atom.mass()*gVector.dot(atom.r);
    }
          
    /**
    * Force exerted by the field on the atom.
    *
    * @return the vector force exerted on the atom
    */
    public Space.Vector force(Atom atom) {
        gForce.Ea1Tv1(atom.mass(),gVector);
        return gForce;
    }
    
    /**
     * main method to test and demonstrate use of this class
     */
    public static void main(String[] args) {
        
        Frame f = new Frame();   //create a window
        f.setSize(600,350);
	    IntegratorVelocityVerlet integrator1 = new IntegratorVelocityVerlet();
	    integrator1.setTimeStep(0.02);  //20 fs time step
	    SpeciesDisks speciesDisks1 = new SpeciesDisks();
	    Phase phase1 = new Phase();
	    P2LennardJones p2LennardJones= new P2LennardJones();  
	    Controller controller1 = new Controller();
		Controller.Button button = controller1.new Button();
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    DisplayBox displayEnergy = new DisplayBox.Energy(phase1); //shows energy to follow conservation (not conserved because of field and PBC)
	    phase1.energy.setActive(true);
	    //unique part for this class
	    phase1.addField(new PotentialFieldGravity(phase1,1e12)); //a huge multiplier is needed to have any noticable effect!

		Simulation.instance.setBackground(Color.yellow);
		Simulation.instance.elementCoordinator.go(); 
        f.add(Simulation.instance);
        f.pack();
        f.show();
        f.addWindowListener(new WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(WindowEvent e) {System.exit(0);}
        });
    }//end of main
       
}



