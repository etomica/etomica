//This class includes a main method to demonstrate its use
package etomica;
import etomica.units.*;

/* History
 * 08/12/03 (DAK) use sim instead of space in AtomFactoryMono constructor
 */
public class SpeciesWalls extends Species implements EtomicaElement {

    private double mass = Default.ATOM_MASS;
    private boolean stationary = false;
/** 
 *  Wall type array.  Each atom has its own type, which specifies its length and orientation.
 *  Examples of use:
 *  Use one molecule, with atoms of different type to set up box, cylinder, etc.
 *  Use several molecules each with one atom to set up parallel walls.
 */
    public AtomType.Wall[] protoType;

    private static AtomFactoryHetero makeFactory(Simulation sim, int nA) {
        AtomFactoryMono[] f = new AtomFactoryMono[nA];
        for(int i=0; i<nA; i++) {
            f[i] = new AtomFactoryMono(sim, sim.iteratorFactory.neighborSequencerFactory());
            AtomType type = new AtomType.Wall(f[i], Default.ATOM_MASS, Double.MAX_VALUE, 0, 0, 0);// arguments are mass, color, length, angle(degrees)  
            f[i].setType(type);
        }
        AtomFactoryHetero fm = new AtomFactoryHetero(sim, f);
        return fm;
    }
        
    /**
    * Default constructor.  Creates species containing 1 molecule with 1 horizontal wall atom.
    */
    public SpeciesWalls() {
        this(Simulation.instance);
    }
    public SpeciesWalls(int n) {
        this(Simulation.instance, n);
    }
    public SpeciesWalls(Simulation sim) {
        this(sim, 1);
    }
    public SpeciesWalls(Simulation sim, int n) {
        this(sim, n, 1);
    }
    public SpeciesWalls(int nM, int nA) {
        this(Simulation.instance, nM, nA);
    }

    public SpeciesWalls(Simulation sim, int nM, int nA) {
        this(sim, nM, nA, Double.MAX_VALUE, 0, 0, 0);
    }
    
    public SpeciesWalls(Simulation sim, int nM, int nA, double length, int xAngle, int yAngle, int zAngle) {  //angle is in radians
        super(sim, makeFactory(sim, nA));
        factory.setSpecies(this);
        //get a handle to the factories used to make the walls and create an array (protoType)
        //of handles to the types attached to each
        AtomFactoryMono[] subfactory = ((AtomFactoryMono[])((AtomFactoryHetero)factory).childFactory());
        protoType = new AtomType.Wall[nA];
        nMolecules = nM;
        for(int i=0; i<nA; i++) {
           AtomType.Wall type = (AtomType.Wall)subfactory[i].type();
//           type.setStationary(true);
           type.setLength(length);
           type.setXAngle(xAngle);
           type.setYAngle(yAngle);
           type.setZAngle(zAngle);
           protoType[i] = type;
        }
        factory.setConfiguration(new SpeciesWalls.ConfigurationParallel(sim));
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species with 'molecules' constructed from wall-shaped 'atoms'");
        return info;
    }

    public void setAlignment(Constants.Alignment a) {
        for(int i=0; i<protoType.length; i++) protoType[i].setAlignment(a);
    }
    public Constants.Alignment getAlignment() {return protoType[0].getAlignment();}
    
    public int getThickness() {return protoType[0].getThickness();}
    public void setThickness(int t) {
        for(int i=0; i<protoType.length; i++) protoType[i].setThickness(t);
    }                    
    
    public double getMass() {return mass;}
    public void setMass(double m) {
        mass = m;
        allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.coord.setMass(mass);}});
    }
    public boolean isStationary() {return stationary;}
    public void setStationary(boolean b) {
        stationary = b;
        allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.coord.setStationary(stationary);}});
    }
    public Dimension getMassDimension() {return Dimension.MASS;}
                
    public double getLength() {return protoType[0].getLength();}
    public void setLength(double d) {
        for(int i=0; i<protoType.length; i++) protoType[i].setLength(d);
    }                    
    public Dimension getLengthDimension() {return Dimension.LENGTH;}
    
    public double getXAngle() {return(protoType[0].getXAngle());}
    public void setXAngle(double t) { protoType[0].setXAngle(t);}
    public Dimension getXAngleDimension() {return Dimension.ANGLE;}
    
    public double getYAngle() {return(protoType[0].getYAngle());}
    public void setYAngle(double t) { protoType[0].setYAngle(t);}
    public Dimension getYAngleDimension() {return Dimension.ANGLE;}
    
    public double getZAngle() {return(protoType[0].getZAngle());}
    public void setZAngle(double t) { protoType[0].setZAngle(t);}
    public Dimension getZAngleDimension() {return Dimension.ANGLE;}
                
    public boolean getLongWall() {return(protoType[0].isLongWall());}
    public void setLongWall(boolean b) { protoType[0].setLongWall(b);}
            
    public double getTemperature() {return(protoType[0].getTemperature());}
    public void setTemperature(double t) {protoType[0].setTemperature(t);}
    public Dimension getTemperatureDimension() {return Dimension.TEMPERATURE;}
            
    public boolean getAdiabatic() { return(protoType[0].isAdiabatic());}
    public void setAdiabatic(boolean a) { protoType[0].setAdiabatic(a);}

    //Class for arranging walls in parallel
    public class ConfigurationParallel extends Configuration {
    
        private double angle;
        private boolean horizontal, vertical;
        private boolean longWall;  //If true, specifies that the wall extends the whole length of the simulation volume
        private double temperature = Default.TEMPERATURE;
        private double placement;
        
        public ConfigurationParallel(Simulation sim) {
            super(sim);
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
        }
        
        public final boolean isLongWall() {return longWall;}
        public final void setLongWall(boolean s) {longWall = s;}
        
        /**
        * Placement of first wall, as a fraction of the distance from the origin to the end of the phase
        */
        public final double getPlacement() {return placement;}
        public final void setPlacement(double p) {placement = p;}

        public void initializePositions(AtomIterator[] iter) {/*no implementation*/}
        /**
        * Sets wall coordinates 
        */
        public void initializeCoordinates(Atom atom) {  //doesn't handle wall that is not either horizontal or vertical
            AtomTreeNodeGroup mnode = (AtomTreeNodeGroup)atom.node;
            AtomIterator iterator = new AtomIteratorList(((AtomTreeNodeGroup)atom.node).childList);
       //     Space.Vector d = m.parentPhase().dimensions();  //what if parentPhase is null?
            double x, y;
            double h = Default.BOX_SIZE;//d.component(1);
            double w = Default.BOX_SIZE;//d.component(0);
            
            //assume long wall
            if(horizontal) {
                x = Double.MIN_VALUE;
                y = placement*Default.BOX_SIZE;//d.component(1);
                w = 100.;//Double.MAX_VALUE;   //crashes with max-value
            }
            else {//vertical
                x = placement*Default.BOX_SIZE;//d.component(0);
                y = Double.MIN_VALUE;
                h = 100.;//Double.MAX_VALUE;  //crashes with max-value
            }
            
            int i = 0;
            double delta;
            double xyNext;
            double wh;
            if(horizontal) {
                delta = (h-y)/(mnode.childAtomCount()-1);
                i = 1;
                xyNext = y;
                wh = w;
            }
            else { //vertical
                delta = (w-x)/(mnode.childAtomCount()-1);
                i = 0;
                xyNext = x;
                wh = h;
            }                    //2D explicit
            iterator.reset();
            while(iterator.hasNext()) {//equally space all "wall atoms"
                Atom a = iterator.next();
                Space.Vector r = a.coord.position();
                a.coord.momentum().E(0.0);
                a.coord.setStationary(stationary);
                r.setX(i,xyNext);
                xyNext += delta;
                ((AtomType.Wall)a.type).setLength(wh);   //length of wall
                ((AtomType.Wall)a.type).setXAngle(angle);
                ((AtomType.Wall)a.type).setTemperature(temperature);
            }
        }//end of initializeCoordinates
                
  }//end of ConfigurationParallel

  /**
   * Method to demonstrate the use of this class.
   * Simulates a system of hard spheres and a wall.
   */
/*    public static void main(String[] args) {

        Simulation sim = new Simulation();
        Simulation.instance = sim;
//	    SpeciesSpheresMono species = new SpeciesSpheresMono();
        SpeciesSpheres species = new SpeciesSpheres();
        species.setNMolecules(5);
        SpeciesWalls walls = new SpeciesWalls();
        walls.setStationary(true);
	    Phase phase = new Phase();
	    IntegratorHard integrator = new IntegratorHard();
	    P2HardSphere potential = new P2HardSphere();
        P2HardSphereWall wallPotential = new P2HardSphereWall(); //hard wall-sphere interaction
	    Controller controller = new Controller();
	    DisplayPhase display = new DisplayPhase();
	    IntegratorMD.Timer timer = integrator.new Timer(integrator.chronoMeter());
	    timer.setUpdateInterval(10);
	//	elementCoordinator.go();

    //    etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        
        //make the wall vertical
        ((ConfigurationParallel)walls.factory.getConfiguration()).setAngle((int)Degree.UNIT.toSim(90.));
        //wall (2nd species added) is set to species 1 automatically, spheres are species 0; 
        //set new potential to be for 0-1 interaction       
 
        Simulation.instance.elementCoordinator.go();
        
        potential.setSpecies(species, species);
        wallPotential.setSpecies(species, walls);
        
        Simulation.makeAndDisplayFrame(Simulation.instance);
        
    }//end of main
    */
}
