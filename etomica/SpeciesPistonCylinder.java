//This class includes a main method to demonstrate its use
package etomica;
import etomica.units.*;
import etomica.units.Dimension;
//Java2 imports
//import java.util.HashMap;

import etomica.utility.HashMap;


import java.awt.Graphics;

/** Four walls arranged as a piston-cylinder apparatus.  All but one
 *  wall (the piston) is stationary.
 *  Defines a Boundary inner class that sets up a phase boundary to coincide
 *  with the dimensions of the apparatus.  To implement, add an instance of the boundary to the phase 
 *  (see main method for an example).
 */
 
public class SpeciesPistonCylinder extends SpeciesWalls implements Space.Boundary.Maker, EtomicaElement {

    public final String getVersion() {return "SpeciesPistonCylinder:01.06.05/"+super.getVersion();}
 
    private Constants.Direction direction = Constants.NORTH;
    Space2D.Vector dimensions = new Space2D.Vector();
    private Atom piston;
    
   /**
    * Vector giving the unit normal to the piston, away from the inside of the cylinder.
    */
    private Space2D.Vector unitNormal = new Space2D.Vector();
    private int thickness = 4;  //cylinder walls thickness in pixels
    private int pistonThickness = 8;//piston thickness
    private double diameter = 28.0;//diameter of cylinder
    private double length = 20.0;//length of cylinder
    private HashMap deployedFields = new HashMap(5);
    
    public SpeciesPistonCylinder() {
        this(Simulation.instance);
    }
    public SpeciesPistonCylinder(Simulation sim) {
        super(sim,1,4);  //1 molecule, 4 atoms
        java.awt.Color color = Constants.DARK_RED;
        double longLength = length;
        double shortLength = diameter;
        double angle = 0;
        protoType[0] = new AtomType.Wall(Default.ATOM_MASS, color, shortLength, angle); //piston
        protoType[1] = new AtomType.Wall(Default.ATOM_MASS, color, longLength, angle); 
        protoType[2] = new AtomType.Wall(Default.ATOM_MASS, color, shortLength, angle); 
        protoType[3] = new AtomType.Wall(Default.ATOM_MASS, color, longLength, angle); 
        protoType[0].setThickness(pistonThickness);
        for(int i=1; i<4; i++) protoType[i].setThickness(thickness);
        protoType[0].setLongWall(false);
        protoType[1].setLongWall(true);
        protoType[2].setLongWall(false);
        protoType[3].setLongWall(true);
        setUnitNormal();
        moleculeConfiguration = new PistonCylinderConfiguration(this);
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
     * Overrides superclass so that "molecule" (piston-cylinder system) can have only 4 "atoms" (3 walls + piston).
     */
    public void setAtomsPerMolecule(int i) {super.setAtomsPerMolecule(4);} //likewise

    /**
     * Makes molecule and adds to phase a PotentialField that applies a pressure to the piston
     */
    protected Molecule makeMolecule(Phase phase) {
        Molecule m = new Molecule(this, phase, protoType);
        PistonPressureField field = (PistonPressureField)deployedFields.get(phase);
        if(field != null) phase.removeField(field);
        field = new PistonPressureField(phase, new Atom.Iterator.Singlet(m.firstAtom()));
        phase.addField(field);
        deployedFields.put(phase,field);
        return m;
    } 
              
  /**
   * Accessor method for the width of the piston chamber.
   */
  public void setDiameter(double d) {
     diameter = d;
     protoType[0].setLength(diameter);
     protoType[2].setLength(diameter);
     moleculeConfiguration.initializeCoordinates();
  }
  /**
   * Accessor method for the width of the piston chamber.
   */
  public double getDiameter() {return diameter;}
  public Dimension getDiameterDimension() {return Dimension.LENGTH;}
  
  /**
   * Accessor method for the length of the piston chamber.
   */
  public void setLength(double l) {
     length = l;
     protoType[1].setLength(length);
     protoType[3].setLength(length);
     moleculeConfiguration.initializeCoordinates();
  }
  /**
   * Accessor method for the length of the piston chamber.
   */
  public double getLength() {return length;}
  public Dimension getLengthDimension() {return Dimension.LENGTH;}
  
  /**
   * Sets a vector giving the unit normal to the piston, away from the inside of the cylinder.
   */
  private final void setUnitNormal() {
    if(direction == Constants.NORTH)      {unitNormal.x =  0.0; unitNormal.y = +1.0;}//N
    else if(direction == Constants.EAST)  {unitNormal.x = -1.0; unitNormal.y =  0.0;}//E
    else if(direction == Constants.SOUTH) {unitNormal.x =  0.0; unitNormal.y = -1.0;}//S
    else                                  {unitNormal.x = +1.0; unitNormal.y =  0.0;}//W
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
    moleculeConfiguration.initializeCoordinates();
  }
  
    
  /**
   *  Accessor method for thickness of the cylinder walls (in pixels) as drawn to the screen.
   *  This value has no effect on the simulation, only its on-screen rendering.
   */
  public final int getThickness() {return thickness;}
  /**
   *  Accessor method for thickness of the cylinder walls (in pixels) as drawn to the screen.
   *  This value has no effect on the simulation, only its on-screen rendering.
   */
  public final void setThickness(int t) {
    thickness = t;
    for(int i=1; i<atomsPerMolecule; i++) {protoType[i].setThickness(t);}
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
  public final class PistonPressureField extends PotentialField implements PotentialField.Soft {
    private double pressure = 0.0;
    private Space.Vector forceVector = parentSimulation().space().makeVector();
    private double force = 0.0;
    public PistonPressureField(Phase p, Atom.Iterator iter) {
        super(p,iter); 
        setPressure(Bar.UNIT.toSim(500.));
    }
    /**
     * Accessor method for the pressure applied to the piston
     */
    public void setPressure(double p) {
        pressure = p;
        force = pressure*diameter; //compute force as pressure * length
        if(phase() != null && phase().integrator() != null && phase().integrator().isInitialized()) 
                    phase().integrator().reset();
    }
    /**
     * Accessor method for the (3D) pressure applied to the piston
     */
    public double getPressure() {return pressure;}
    /**
     * Returns dimensions of (3D) pressure
     */
    public Dimension getPressureDimension() {return Dimension.PRESSURE;}
    /**
     * Force as computed from current pressure and piston diameter, and direction of piston-cylinder system.
     */
    public Space.Vector force(Atom a) {
        forceVector.Ea1Tv1(force,unitNormal);
        return forceVector;
    }
    /**
     * Always returns zero
     */
    public double energy(Atom a) {return 0.0;}
    
  }//end of PistonPressureField
  
    
  /**
   * Molecule.Configuration class that sets up all the walls of the piston-cylinder.
   * Accounts for the direction (N/E/S/W) of the apparatus.
   */
  private final class PistonCylinderConfiguration extends Molecule.Configuration {
      
      PistonCylinderConfiguration(SpeciesPistonCylinder s) {
        super(s);
      }
      public void initializeCoordinates(Molecule m) {
        m.atomIterator.reset();
        Atom[] atoms = new Atom[4];
        double wallThickness = thickness/BaseUnit.Length.Sim.TO_PIXELS;
        double pistonWallThickness = pistonThickness/BaseUnit.Length.Sim.TO_PIXELS;
        int i=0;
        while(m.atomIterator.hasNext()) {atoms[i++] = m.atomIterator.next();}
        atoms[0].setStationary(false); //piston
 //       atoms[0].setStationary(true); //piston
        atoms[1].setStationary(true);
        atoms[2].setStationary(true);
        atoms[3].setStationary(true);
        if(((SpeciesPistonCylinder)parentSpecies).direction == Constants.NORTH) {
            protoType[0].setAlignment(Constants.HORIZONTAL);
            protoType[1].setAlignment(Constants.VERTICAL);
            protoType[2].setAlignment(Constants.HORIZONTAL);
            protoType[3].setAlignment(Constants.VERTICAL);
            atoms[0].r.setComponent(0,0.0); //top (piston)
            atoms[0].r.setComponent(1,0.0);
            atoms[1].r.setComponent(0,diameter-wallThickness); //right wall
            atoms[1].r.setComponent(1,0.0);
            atoms[2].r.setComponent(0,0.0); //bottom
            atoms[2].r.setComponent(1,length-wallThickness);
            atoms[3].r.E(0.0); //left wall
        }
        else if(((SpeciesPistonCylinder)parentSpecies).direction == Constants.EAST) {
            protoType[3].setAlignment(Constants.HORIZONTAL);
            protoType[0].setAlignment(Constants.VERTICAL);
            protoType[1].setAlignment(Constants.HORIZONTAL);
            protoType[2].setAlignment(Constants.VERTICAL);
            atoms[3].r.E(0.0);
            atoms[0].r.setComponent(0,length-pistonWallThickness);
            atoms[0].r.setComponent(1,0.0);
            atoms[1].r.setComponent(0,0.0);
            atoms[1].r.setComponent(1,diameter-wallThickness);
            atoms[2].r.E(0.0);
        }
        else if(((SpeciesPistonCylinder)parentSpecies).direction == Constants.SOUTH) {
            protoType[2].setAlignment(Constants.HORIZONTAL);
            protoType[3].setAlignment(Constants.VERTICAL);
            protoType[0].setAlignment(Constants.HORIZONTAL);
            protoType[1].setAlignment(Constants.VERTICAL);
            atoms[2].r.E(0.0);
            atoms[3].r.setComponent(0,diameter-wallThickness);
            atoms[3].r.setComponent(1,0.0);
            atoms[0].r.setComponent(0,0.0);
            atoms[0].r.setComponent(1,length-pistonWallThickness);
            atoms[1].r.E(0.0);
        }
        else { //WEST
            protoType[1].setAlignment(Constants.HORIZONTAL);
            protoType[2].setAlignment(Constants.VERTICAL);
            protoType[3].setAlignment(Constants.HORIZONTAL);
            protoType[0].setAlignment(Constants.VERTICAL);
            atoms[1].r.E(0.0);
            atoms[2].r.setComponent(0,length-wallThickness);
            atoms[2].r.setComponent(1,0.0);
            atoms[3].r.setComponent(0,0.0);
            atoms[3].r.setComponent(1,diameter-wallThickness);
            atoms[0].r.E(0.0);
        }
        computeDimensions();
      }//end of initializeCoordinates
      
      public void computeDimensions() {SpeciesPistonCylinder.this.computeDimensions();}
      
  }//end of PistonCylinderConfiguration
  
        public void computeDimensions() {
            double d = diameter - thickness / BaseUnit.Length.Sim.TO_PIXELS;
            double l = length - thickness / BaseUnit.Length.Sim.TO_PIXELS;
            dimensions.x = (direction == Constants.NORTH || direction == Constants.SOUTH) ? d : l;
            dimensions.y = (direction == Constants.NORTH || direction == Constants.SOUTH) ? l : d;
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
            Molecule m = getAgent(this.phase()).firstMolecule();
            int index20 = (direction==Constants.NORTH || direction==Constants.SOUTH) ? 1 : 0;
            int index31 = 1-index20;//0 or 1, opposite of index20
              //volume is diameter of cylinder times distance between piston (first atom) and base (atom 2) of cylinder
 //           return (diameter - thickness / BaseUnit.Length.Sim.TO_PIXELS)
            return Math.abs(m.getAtom(3).r.component(index31) - m.getAtom(1).r.component(index31))
                    * Math.abs(m.firstAtom().r.component(index20) - m.getAtom(2).r.component(index20));
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
            Species.Agent agent = SpeciesPistonCylinder.this.getAgent(phase);
            piston = agent.firstAtom();
            setLabel("Hold/Release piston");
            this.phase = phase;
        }
        
        public void actionPerformed() {
            piston.setStationary(!piston.isStationary()); //set stationary
            piston.p.E(0.0);                              //set momentum to zero
            phase.integrator().reset();                   //update integrator
        }
    }
    
    /**
     * Demonstrates and tests this class.
     */
    public static void main(String[] args) {
        javax.swing.JFrame f = new javax.swing.JFrame();
        f.setSize(600,350);
	    SpeciesDisks speciesDisks1 = new SpeciesDisks();
	    Phase phase1 = new Phase();
	    P2HardDisk P2HardDisk1 = new P2HardDisk();
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
		Simulation.instance.setBackground(java.awt.Color.yellow);
        
        //part unique to this class
	    IntegratorHard integratorHard1 = new IntegratorHardField();
        SpeciesPistonCylinder pistonCylinder = new SpeciesPistonCylinder();
        pistonCylinder.setLength(20.);
        ((DisplayPhase)Simulation.instance.display(0)).setColorScheme(new ColorSchemeByType());
        P2DiskWall p2DiskWall = new P2DiskWall();
        p2DiskWall.setSpeciesIndex(0,1);
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
        DeviceSlider pressureSlider = new DeviceSlider((PistonPressureField)phase1.firstField(),"pressure");
        Simulation.instance.add(pressureSlider.graphic(null));
        pressureSlider.setUnit(new Unit(Bar.UNIT));
        pressureSlider.setMinimum(50);
        pressureSlider.setMaximum(1000);
        //end of unique part

		f.getContentPane().add(Simulation.instance);
		f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
  
}
