//This class includes a main method to demonstrate its use
package etomica;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.Iterator;
import etomica.units.*;

public class SpeciesWalls extends Species implements EtomicaElement {

    public static String version() {return "01.03.05.0";}

/** 
 *  Wall type array.  Each atom has its own type, which specifies its length and orientation.
 *  Examples of use:
 *  Use one molecule, with atoms of different type to set up box, cylinder, etc.
 *  Use several molecules each with one atom to set up parallel walls.
 */
    public AtomType.Wall[] protoType;

    /**
    * Default constructor.  Creates species containing 1 molecule with 1 horizontal wall atom.
    */
    public SpeciesWalls() {
        this(Simulation.instance);
    }
    public SpeciesWalls(Simulation sim) {
        this(sim,1,1,Double.MAX_VALUE,0);
    }
    
    /**
     * Sets up nM sets of wall 'molecules' with each molecule having nA wall 'atoms'.
     * Defines but does not fill protoType array.
     */
    public SpeciesWalls(int nM, int nA) {
        this(Simulation.instance);
    }
    public SpeciesWalls(Simulation sim, int nM, int nA) {
        super(sim);
        protoType = new AtomType.Wall[nA];
        atomsPerMolecule = nA;
        nMolecules = nM;
    }
    
    public SpeciesWalls(Simulation sim, int nM, int nA, double length, int angle) {  //angle is in radians
        this(sim, nM, nA);
        for(int i=0; i<nA; i++) {protoType[i] = new AtomType.Wall(Default.ATOM_MASS, Default.ATOM_COLOR, length, angle);}  // arguments are mass, color, length, angle(degrees)
        setStationary(true);

        moleculeConfiguration = new SpeciesWalls.ConfigurationParallel(this);
    }

    public SpeciesWalls(int nM, AtomType.Wall[] type) {
        this(Simulation.instance, nM, type);
    }
    public SpeciesWalls(Simulation sim, int nM, AtomType.Wall[] type) { 
        super(sim);
        protoType = type;
        atomsPerMolecule = type.length;
        nMolecules = nM;
        setStationary(true);

        moleculeConfiguration = new SpeciesWalls.ConfigurationParallel(this);
        moleculeConfiguration.setParentSpecies(this);
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species with 'molecules' constructed from wall-shaped 'atoms'");
        return info;
    }

    protected Molecule makeMolecule(Phase phase) {
        return new Molecule(this, phase, protoType);
    } 
              
        public final void setAlignment(Constants.Alignment a) {
            protoType[0].setAlignment(a);
        }
        public final Constants.Alignment getAlignment() {return protoType[0].getAlignment();}
    // Exposed Properties --- not implemented now because they must tie to AtomType array
    
/*  public final int getThickness() {return ((AtomWall)firstAtom()).getThickness();}
  public final void setThickness(int t) {((AtomWall)firstAtom()).setThickness(t);}
    
    public final double getMass() {return protoType.mass();}
    public final void setMass(double mass) {protoType.setMass(mass);}
                
    public final double getLength() {return protoType.length();}
    public void setLength(double d) {protoType.setLength(d);}
                    
    public final Color getColor() {return protoType.color();}
    public final void setColor(Color c) {protoType.setColor(c);}*/
    
    //Class for arranging walls in parallel
    public class ConfigurationParallel extends Molecule.Configuration {
    
        private double angle;
        private boolean horizontal, vertical;
        private boolean longWall;  //If true, specifies that the wall extends the whole length of the simulation volume
        private double temperature = Default.TEMPERATURE;
        private double placement;
        
        public ConfigurationParallel(Species parent) {
            super(parent);
            parentSpecies = parent;
            setAngle(0.0);
            setLongWall(false);
            setPlacement(0.0);
        }
          
        public final double getAngle() {return angle;}
        public final void setAngle(double t) {
            double pi = Math.PI;
            t = (Math.abs(t) > pi/4.) ? pi/2. : 0.0;  //For now, allow only values for vertical or horizontal walls
            angle = (t <= 2.0*pi) ? t : (t % (2.0*pi));
            horizontal = (angle == 0.0) || (Math.abs(angle) == pi);
            vertical = (Math.abs(angle) == pi/2.) || (Math.abs(angle) == 1.5*pi);
            initializeCoordinates();
        }
        
        public final double getTemperature() {return temperature;}
        public final void setTemperature(double t) {
            temperature = t;
            initializeCoordinates();
        }
        
        public final boolean isLongWall() {return longWall;}
        public final void setLongWall(boolean s) {
            longWall = s;
            initializeCoordinates();
        }
        
        /**
        * Placement of first wall, as a fraction of the distance from the origin to the end of the phase
        */
        public final double getPlacement() {return placement;}
        public final void setPlacement(double p) {placement = p;}

        /**
        * Sets wall coordinates 
        */
        public void initializeCoordinates(Molecule m) {  //doesn't handle wall that is not either horizontal or vertical
            Space.Vector d = m.parentPhase().dimensions();  //what if parentPhase is null?
            double x, y;
            double h = d.component(1);
            double w = d.component(0);
            
            //assume long wall
            if(horizontal) {
                x = Double.MIN_VALUE;
                y = placement*d.component(1);
                w = 100.;//Double.MAX_VALUE;   //crashes with max-value
            }
            else {//vertical
                x = placement*d.component(0);
                y = Double.MIN_VALUE;
                h = 100.;//Double.MAX_VALUE;  //crashes with max-value
            }
            
            int i = 0;
            double delta;
            double xyNext;
            double wh;
            if(horizontal) {
                delta = (h-y)/(m.atomCount-1);
                i = 1;
                xyNext = y;
                wh = w;
            }
            else { //vertical
                delta = (w-x)/(m.atomCount-1);
                i = 0;
                xyNext = x;
                wh = h;
            }                    //2D explicit
            m.atomIterator.reset();
            while(m.atomIterator.hasNext()) {//equally space all "wall atoms"
                Atom a = m.atomIterator.next();
                Space.Vector r = a.coordinate.position();
                a.coordinate.momentum().E(0.0);
                r.setComponent(i,xyNext);
                xyNext += delta;
                ((AtomType.Wall)a.type).setLength(wh);   //length of wall
                ((AtomType.Wall)a.type).setAngle(angle);
                ((AtomType.Wall)a.type).setTemperature(temperature);
            }
        }//end of initializeCoordinates
        
        protected void computeDimensions() {
            if(parentSpecies()==null) return;
    /*       Molecule m = parentSpecies.getMolecule();
    //        initializeCoordinates(m);
            if(horizontal) {
    //            dim[0] = ((AtomType.Wall)m.firstAtom().type).getLength();
                dim[0] = Double.MAX_VALUE;
                dim[1] = 0.0;
            }
            else if(vertical) {
                dim[0] = 0.0;
                dim[1] = Double.MAX_VALUE;
    //            dim[1] = ((AtomType.Wall)m.firstAtom().type).getLength();
            }
            else {
                //does not handle walls that are neither horizontal nor vertical
            }*/
        }//end of computeDimensions
        
  }//end of ConfigurationParallel

  /**
   * Method to demonstrate the use of this class.
   * Simulates a system of hard disks and a wall.
   */
    public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        Simulation.makeSimpleSimulation();  
        
        //here's the part unique to this class
        SpeciesWalls walls = new SpeciesWalls();
        //make the wall vertical
        ((SpeciesWalls.ConfigurationParallel)walls.moleculeConfiguration).setAngle(Degree.UNIT.toSim(90.));
        P2DiskWall wallPotential = new P2DiskWall(); //hard wall-disk interaction
        //wall (2nd species added) is set to species 1 automatically, disks are species 0; 
        //set new potential to be for 0-1 interaction       
        wallPotential.setSpecies2Index(1);  //could set either index (species1Index or species2Index); both are 0 by default
        //end of unique part
 
        Simulation.instance.elementCoordinator.go();
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
}
