//This class includes a main method to demonstrate its use
package etomica;
import etomica.units.*;
import etomica.units.Dimension;

/** 
 * Four walls arranged as a piston-cylinder apparatus.  All but one
 * wall (the piston) is stationary.
 * Defines a Boundary inner class that sets up a phase boundary to coincide
 * with the dimensions of the apparatus.  To implement, add an instance of the boundary to the phase 
 * (see main method for an example).
 *
 * @author David Kofke
 */
 
public class SpeciesPistonCylinder extends SpeciesWalls implements Space.Boundary.Maker, EtomicaElement, Constants {

    public String getVersion() {return "SpeciesPistonCylinder:01.06.05/"+super.getVersion();}
 
    //making static as a hack until direction can be associated with each molecule
    public static Constants.Direction direction = TOP;
    Space2D.Vector dimensions = new Space2D.Vector();
    private Atom piston;
    
   /**
    * Vector giving the unit normal to the piston, away from the inside of the cylinder.
    */
    private Space2D.Vector unitNormal = new Space2D.Vector();
    private int thickness = 4;  //cylinder walls thickness in pixels
    private int pistonThickness = 8;//piston thickness
    private double diameter = 30.0;//diameter of cylinder
    private double length = 20.0;//length of cylinder
    
    public SpeciesPistonCylinder() {
        this(Simulation.instance);
    }
    public SpeciesPistonCylinder(Simulation sim) {
        super(sim,1,4);  //1 molecule, 4 atoms
        double longLength = length;
        double shortLength = diameter;
        double angle = 0;
        protoType[0].setLength(shortLength);
        protoType[1].setLength(longLength);
        protoType[2].setLength(shortLength);
        protoType[3].setLength(longLength);
        for(int i=0; i<4; i++) {
            protoType[i].setThickness(thickness);
        }
        protoType[0].setThickness(pistonThickness);
        protoType[0].setLongWall(false);
        protoType[1].setLongWall(true);
        protoType[2].setLongWall(false);
        protoType[3].setLongWall(true);
        setUnitNormal();
        factory.setSpecies(this);
        factory.setConfiguration(new PistonCylinderConfiguration(sim));
        computeDimensions();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Three fixed and one movable wall in a piston-cylinder configuration");
        return info;
    }

    /**
     * Overrides superclass so that only 1 "molecule" can be constructed.
     */
    public void setNMolecules(int i) {super.setNMolecules(1);}  //override so nMolecules cannot be changed
              
  /**
   * Mutator method for the width of the piston chamber.
   */
  public void setDiameter(double d) {
     diameter = d;
     protoType[0].setLength(diameter);
     protoType[2].setLength(diameter);
  }
  /**
   * Accessor method for the width of the piston chamber.
   */
  public double getDiameter() {return diameter;}
  public Dimension getDiameterDimension() {return Dimension.LENGTH;}
  
  /**
   * Mutator method for the length of the piston chamber.
   */
  public void setLength(double l) {
     length = l;
     protoType[1].setLength(length);
     protoType[3].setLength(length);
  }
  /**
   * Accessor method for the length of the piston chamber.
   */
  public double getLength() {return length;}
  public Dimension getLengthDimension() {return Dimension.LENGTH;}
  
  /**
   * Sets a vector giving the unit normal to the piston, toward the inside of the cylinder.
   */
  private final void setUnitNormal() {
    if(direction == TOP)      {unitNormal.setX(0,0.0); unitNormal.setX(1,+1.0);}//N
    else if(direction == RIGHT)  {unitNormal.setX(0,-1.0); unitNormal.setX(1, 0.0);}//E
    else if(direction == BOTTOM) {unitNormal.setX(0, 0.0); unitNormal.setX(1, -1.0);}//S
    else /*LEFT*/                         {unitNormal.setX(0,+1.0); unitNormal.setX(1,0.0);}//W
  }
  
  /**
   * Direction of orientation of the piston-cylinder system.  
   * Input is a compass direction (e.g., Constants.NORTH) that indicates the placement of the piston.
   * For example, NORTH has the piston on top, EAST has it on the right, etc.
   */
  public Constants.Direction getDirection() {return direction;}
  /**
   * Direction of orientation of the piston-cylinder system.  
   * Input is a compass direction (e.g., Constants.NORTH) that indicates the placement of the piston.
   * For example, NORTH has the piston on top, EAST has it on the right, etc.
   * Default is NORTH.
   */
  public void setDirection(Constants.Direction d) {
    direction = d;
    setUnitNormal();
  }
  
    
  /**
   *  Accessor method for thickness of the piston (in pixels) as drawn to the screen.
   *  This value has no effect on the simulation, only its on-screen rendering.
   */
  public final int getPistonThickness() {return pistonThickness;}
  /**
   *  Accessor method for thickness of the piston (in pixels) as drawn to the screen.
   *  This value has no effect on the simulation, only its on-screen rendering.
   */
  public final void setPistonThickness(int t) {
    pistonThickness = t;
    protoType[0].setThickness(t);
  }
  
  //interface method
  /**
   * Returns a new Boundary object suitable for use with this species.
   * Boundary is not periodic.
   */
  public Space.Boundary makeBoundary(Space.Boundary.Type t) {
    if(t != BOUNDARY) {throw new IllegalArgumentException();}
    return this.new Boundary();
  }
  /**
   * Returns a single-element array containing the type indicating a Piston-cylinder boundary.
   */
  public Space.Boundary.Type[] boundaryTypes() {return new Space.Boundary.Type[] {BOUNDARY};}
  /**
   * Returns true.
   */
  public boolean requiresSpecialBoundary() {return true;}
  public static final Space.Boundary.Type BOUNDARY = 
    new Space2D.Boundary.Type("Piston-cylinder") {
        public Constants.TypedConstant[] choices() {return new Constants.TypedConstant[] {this};}
    };
  
  /**
   * Field that applies a constant force against the piston.
   * Direction of force is such that it pushes the piston into the cylinder.
   */
  public final class PistonPressureField extends Potential1Soft implements PotentialHard {
    private double pressure = 0.0;
    private Space.Vector gradientVector = parentSimulation().space().makeVector();
    private double force = 0.0;
    public PistonPressureField() {this(Simulation.instance.hamiltonian.potential);}
    public PistonPressureField(PotentialGroup parent) {
        super(parent); 
        setPressure(Bar.UNIT.toSim(100.));
        iterator = new PistonIterator();
        setSpecies(SpeciesPistonCylinder.this);
    }
    /**
     * Accessor method for the pressure applied to the piston
     */
    public void setPressure(double p) {
        pressure = p;
        force = pressure*diameter; //compute force as pressure * length
//        System.out.println("pressure, diameter, force: "+pressure+" "+diameter+" "+force);
        parentSimulation().resetIntegrators();
//        if(phase() != null && phase().integrator() != null && phase().integrator().isInitialized()) 
//                    phase().integrator().reset();
    }
    /**
     * Accessor method for the pressure applied to the piston
     */
    public double getPressure() {return pressure;}
    /**
     * Returns dimensions of (3D) pressure
     */
    public Dimension getPressureDimension() {return Dimension.PRESSURE;}
    /**
     * Potential gradient as computed from current pressure and piston diameter, 
     * and direction of piston-cylinder system.
     */
    public Space.Vector gradient(Atom a) {
        gradientVector.Ea1Tv1(-force,unitNormal);
 //       System.out.println(gradientVector.toString());
        return gradientVector;
    }
    /**
     * Always returns zero
     */
    public double energy(Atom a) {return 0.0;}
    
    private final class PistonIterator extends AtomIteratorSinglet {
        //Basis will be set to species agent for PistonCylinder.  Assuming
        //there is only one piston-cylinder "molecule", the piston will be
        //the first leaf atom of this agent.  This is the only atom acted
        //on by this potential.
        public void setBasis(Atom a) {
            setAtom(((AtomGroup)a).node.firstLeafAtom());
        }
    }
    
    //implementation of PotentialHard interface, which is done because this field is sometimes
    //used with a hard-potential, constant-field integrator.
    public void bump(AtomPair pair) {}
    
    public double lastCollisionVirial() {return 0.0;}
    
    public Space.Tensor lastCollisionVirialTensor() {return null;}
  }//end of PistonPressureField
  
    
  /**
   * Molecule.Configuration class that sets up all the walls of the piston-cylinder.
   * Accounts for the direction (N/E/S/W) of the apparatus.
   */
  private final class PistonCylinderConfiguration extends Configuration {
      
      PistonCylinderConfiguration(Simulation sim) {
        super(sim);
      }
      public void initializePositions(AtomIterator[] iter) {/*no implementation*/}
      public void initializeCoordinates(Atom atom) {
        AtomIterator iterator = new AtomIteratorTree(atom);
        Atom[] atoms = new Atom[4];
        double wallThickness = thickness/BaseUnit.Length.Sim.TO_PIXELS;
        double pistonWallThickness = pistonThickness/BaseUnit.Length.Sim.TO_PIXELS;
        int i=0;
        iterator.reset();
        while(iterator.hasNext()) {atoms[i++] = iterator.next();}
        atoms[0].coord.setStationary(false); //piston
 //       atoms[0].setStationary(true); //piston
        atoms[1].coord.setStationary(true);
        atoms[2].coord.setStationary(true);
        atoms[3].coord.setStationary(true);
        if(SpeciesPistonCylinder.direction == TOP) {
            protoType[0].setAlignment(HORIZONTAL);
            protoType[1].setAlignment(VERTICAL);
            protoType[2].setAlignment(HORIZONTAL);
            protoType[3].setAlignment(VERTICAL);
            atoms[0].coord.position().setX(0,0.0); //top (piston)
            atoms[0].coord.position().setX(1,0.0);
            atoms[1].coord.position().setX(0,diameter-wallThickness); //right wall
            atoms[1].coord.position().setX(1,0.0);
            atoms[2].coord.position().setX(0,0.0); //bottom
            atoms[2].coord.position().setX(1,length-wallThickness);
            atoms[3].coord.position().E(0.0); //left wall
        }
        else if(SpeciesPistonCylinder.direction == RIGHT) {
            protoType[3].setAlignment(HORIZONTAL);
            protoType[0].setAlignment(VERTICAL);
            protoType[1].setAlignment(HORIZONTAL);
            protoType[2].setAlignment(VERTICAL);
            atoms[3].coord.position().E(0.0);
            atoms[0].coord.position().setX(0,length-pistonWallThickness);
            atoms[0].coord.position().setX(1,0.0);
            atoms[1].coord.position().setX(0,0.0);
            atoms[1].coord.position().setX(1,diameter-wallThickness);
            atoms[2].coord.position().E(0.0);
        }
        else if(SpeciesPistonCylinder.direction == BOTTOM) {
            protoType[2].setAlignment(HORIZONTAL);
            protoType[3].setAlignment(VERTICAL);
            protoType[0].setAlignment(HORIZONTAL);
            protoType[1].setAlignment(VERTICAL);
            atoms[2].coord.position().E(0.0);
            atoms[3].coord.position().setX(0,diameter-wallThickness);
            atoms[3].coord.position().setX(1,0.0);
            atoms[0].coord.position().setX(0,0.0);
            atoms[0].coord.position().setX(1,length-pistonWallThickness);
            atoms[1].coord.position().E(0.0);
        }
        else { //WEST
            protoType[1].setAlignment(HORIZONTAL);
            protoType[2].setAlignment(VERTICAL);
            protoType[3].setAlignment(HORIZONTAL);
            protoType[0].setAlignment(VERTICAL);
            atoms[1].coord.position().E(0.0);
            atoms[2].coord.position().setX(0,length-wallThickness);
            atoms[2].coord.position().setX(1,0.0);
            atoms[3].coord.position().setX(0,0.0);
            atoms[3].coord.position().setX(1,diameter-wallThickness);
            atoms[0].coord.position().E(0.0);
        }
        computeDimensions();
      }//end of initializeCoordinates
      
      public void computeDimensions() {SpeciesPistonCylinder.this.computeDimensions();}
      
  }//end of PistonCylinderConfiguration
  
        public void computeDimensions() {
            double d = diameter - thickness / BaseUnit.Length.Sim.TO_PIXELS;
            double l = length - thickness / BaseUnit.Length.Sim.TO_PIXELS;
            dimensions.setX(0, (direction == TOP || direction == BOTTOM) ? d : l);
            dimensions.setX(1, (direction == TOP || direction == BOTTOM) ? l : d);
        }//end of computeDimensions
      
    /**
     * A phase boundary that coincides with the dimensions of the piston-cylinder system.
     * Dimension of boundary are given by fixed length and diameter of cylinder, while volume
     * of boundary fluctuates according to the position of the piston.
     */
    public class Boundary extends Space2D.BoundaryNone {
        
        public Boundary() {super();}
        public Boundary(Phase p) {super(p);}
        public Space.Vector dimensions() {return SpeciesPistonCylinder.this.dimensions;}
        public double volume() {
            AtomTreeNodeGroup node = (AtomTreeNodeGroup)getAgent(this.phase()).firstMolecule().node;
            int index20 = (direction==TOP || direction==BOTTOM) ? 1 : 0;
            int index31 = 1-index20;//0 or 1, opposite of index20
              //volume is diameter of cylinder times distance between piston (first atom) and base (atom 2) of cylinder
 //           return (diameter - thickness / BaseUnit.Length.Sim.TO_PIXELS)
            double dx = node.getAtom(3).coord.position().x(index31) - node.getAtom(1).coord.position().x(index31);
            double dy = node.getAtom(0).coord.position().x(index20) - node.getAtom(2).coord.position().x(index20);
//            System.out.println(dx + " " + dy);
            return Math.abs(dx*dy);
        }
    }//end of Boundary class
    
    /**
     * An extension of the Action class that toggles the piston between
     * being stationary (fixed) and freely movable.
     */
    public class ActionFixPiston extends etomica.Action {
        Atom piston;
        Phase phase;
        public ActionFixPiston(Phase phase) {
            SpeciesAgent agent = SpeciesPistonCylinder.this.getAgent(phase);
            piston = agent.node.firstLeafAtom();
            setLabel("Hold/Release piston");
            this.phase = phase;
        }
        
        public void actionPerformed() {
            piston.coord.setStationary(!piston.coord.isStationary()); //set stationary
            piston.coord.momentum().E(0.0);                        //set momentum to zero
            phase.integrator().reset();                   //update integrator
        }
    }
    
    /**
     * Demonstrates and tests this class.
     */
/*    public static void main(String[] args) {
        Simulation sim = new Simulation();
        Simulation.instance = sim;
        
	    SpeciesSpheres speciesSpheres1 = new SpeciesSpheres(20);
	    Phase phase1 = new Phase();
	    P2HardSphere P2HardSphere1 = new P2HardSphere();
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
		Simulation.instance.panel().setBackground(java.awt.Color.yellow);
        
        //part unique to this class
	    IntegratorHard integratorHard1 = new IntegratorHardField();
	  //  integratorHard1.setTimeStep(0.01);
        SpeciesPistonCylinder pistonCylinder = new SpeciesPistonCylinder();
        pistonCylinder.setLength(20.);
        ((DisplayPhase)Simulation.instance.display(0)).setColorScheme(new ColorSchemeByType());
        P2HardSphereWall wallPotential = new P2HardSphereWall();
        
        PistonPressureField pressureField = pistonCylinder.new PistonPressureField();

        DeviceSlider pressureSlider = new DeviceSlider(pressureField,"pressure");
        pressureSlider.setUnit(new Unit(Bar.UNIT));
        pressureSlider.setMinimum(50);
        pressureSlider.setMaximum(1000);
        
        Simulation.instance.phase(0).setBoundary(pistonCylinder.new Boundary(Simulation.instance.phase(0))); //have piston-cylinder system define boundary of phase
        Meter thermometer = new MeterTemperature();
        thermometer.setPhase(phase1);
        DisplayBox tBox = new DisplayBox();
        tBox.setMeter(thermometer);
        tBox.setUnit(new Unit(Kelvin.UNIT));
        displayPhase1.setAlign(1,DisplayPhase.BOTTOM);
 
        
	    IntegratorMD.Timer timer = integratorHard1.new Timer(integratorHard1.chronoMeter());
	    timer.setUpdateInterval(10);
		Simulation.instance.elementCoordinator.go();

        DeviceButton fixPiston = new DeviceButton(pistonCylinder.new ActionFixPiston(phase1));

  //    wall-sphere potential
    //order of specification of atom iterators does matter!
        wallPotential.setSpecies(speciesSpheres1, pistonCylinder);
        P2HardSphere1.setSpecies(speciesSpheres1, speciesSpheres1);
        //need to add pressure field iterator
        
    /*    Potential2.Agent potentialAgent = (Potential2.Agent)wallPotential.getAgent(phase1);
        potentialAgent.setIterator(new AtomPairIterator(phase1,
                speciesSpheres1.makeAtomIterator(phase1), pistonCylinder.makeAtomIterator(phase1)));

        //sphere-sphere potential
        potentialAgent = (Potential2.Agent)P2HardSphere1.getAgent(phase1);
        potentialAgent.setIterator(new AtomPairIterator(phase1,
                speciesSpheres1.makeAtomIterator(phase1), speciesSpheres1.makeAtomIterator(phase1)));

        //pressure field on piston
        Potential1.Agent potential1Agent = (Potential1.Agent)pressureField.getAgent(phase1);
        potential1Agent.setIterator(new AtomIteratorSinglet(((SpeciesAgent)pistonCylinder.getAgent(phase1)).firstLeafAtom()));
   * /     
		Simulation.instance.elementCoordinator.go();
//		pistonCylinder.getAgent(phase1).firstLeafAtom().coord.momentum().E(0.0);
		pistonCylinder.setStationary(true);
		pistonCylinder.getAgent(phase1).firstLeafAtom().coord.setStationary(false);
        Simulation.makeAndDisplayFrame(Simulation.instance);
        
  /* //code that displays collision times with atoms 
        final Phase phase = phase1;
	    displayPhase1.addDrawable(new DisplayPhase.Drawable() {
	        public void draw(java.awt.Graphics g, int[] origin, double s) {
                double toPixels = etomica.units.BaseUnit.Length.Sim.TO_PIXELS*s;
                int i=0;
	            for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {
	                IntegratorHardAbstract.Agent agent = (IntegratorHardAbstract.Agent)a.ia;
	                if(agent == null) return;
	                String text = Float.toString((float)agent.collisionTime);
                    Space.Vector r = a.coord.position();
                    int xP = origin[0] + (int)(toPixels*(r.x(0)));
                    int yP = origin[1] + (int)(toPixels*(r.x(1)));
                    g.setColor(java.awt.Color.gray);
	                g.drawString(text, xP, yP-20);
	                g.setColor(java.awt.Color.red);
	                g.drawString(Integer.toString(a.index()), xP-20, yP-20);
	            }
	        }
	    });  * /
    }//end of main
  */
}
