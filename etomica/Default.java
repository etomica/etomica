package etomica;

import etomica.units.Bar;
import etomica.units.Dimension;
import etomica.units.Kelvin;
import etomica.units.UnitSystem;

/**
 * Class holding fields that define many of the default values used in building
 * a simulation.
 * @author kofke
 */

/* History
 * 09/02/03 (DAK) added DO_SLEEP, used by Integrator
 */

public final class Default {
    
 //   public static String WORKING_DIRECTORY = getWorkingDirectory();
 //   public static String WORKING_DIRECTORY = "D:\\Etomica\\";
   public static String WORKING_DIRECTORY = "";
    
    public static String CLASS_DIRECTORY = WORKING_DIRECTORY + "etomica";
    
    public static double ATOM_SIZE = 3.0;  //Angstroms
    
    public static double ATOM_MASS = 40.0; //Daltons
    
    public static int MOLECULE_COUNT = 20;
    
    public static double BOX_SIZE = 30.0;  //Angstroms
    
    public static double TEMPERATURE = Kelvin.UNIT.toSim(300.);
    
    public static double PRESSURE = Bar.UNIT.toSim(1.0);
    
    public static double POTENTIAL_WELL = Kelvin.UNIT.toSim(300.);
    
    public static double POTENTIAL_CUTOFF_FACTOR = 2.5; //dimensionless multiplier for cutoff distance of potential
    
    public static double TIME_STEP = 0.05;  //picoseconds 
    
    public static int HISTORY_PERIOD = 100;
    
    public static boolean TRUNCATE_POTENTIALS = true;
    
    public static boolean IS_GRAPHIC = false;
    
    public static boolean FIX_OVERLAP = false;
    
    public static boolean AUTO_REGISTER = true;
    
    /**
     * Flag indicating how iteration is done in potential calculations.  If
     * true, then hasNext()/next() construct is used, if false, all() construct
     * is used.
     */
    public static final boolean EXPLICIT_LOOP = true;
    
    /**
     * Default value for doSleep field in Integrator class.  The default defined
     * here is <false>, indicating that Integrator should not pause during
     * integration loop.  Instantiation of any SimulationGraphic class changes
     * this default to true, which is appropriate for simulations that use
     * interactive graphics -- doSleep <true> is needed for a responsive
     * interface.  Change in this default has no effect on any Integrators
     * previously constructed.
     */
    public static boolean DO_SLEEP = false;
    
    /**
     * Integer array indicating the maximum number of atoms at each depth in the
     * atom hierarchy.  Maximum depth is given by the size of the array.  Each
     * element of array is log2 of the maximum number of child atoms permitted
     * under one atom.  This is used to assign index values to each atom when it
     * is made.  The indexes permit quick comparison of the relative ordering
     * and/or hierarchical relation of any two atoms.  Sum of digits in array
     * should not exceed 31. For example, {5, 16, 7, 3} indicates 31
     * speciesAgents maximum, 65,536 molecules each, 128 groups per molecule, 8
     * atoms per group (number of species agents is one fewer, because index 0
     * is assigned to species master)
     */
    // powers of 2, for reference:
    //  n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |  14  |  15  |  16  |  17   |  18   |  19   |   20    |
    // 2^n| 2 | 4 | 8 | 16| 32| 64|128|256|512|1024|2048|4096|8192|16,384|32,768|65,536|131,072|262,144|524,288|1,048,576|
    // {speciesRoot, phases, species, molecules, groups, atoms}
    public static int[] BIT_LENGTH = new int[] {1, 4, 4, 14, 6, 2};

    
    /**
     * Sets default atom size, mass, and potential-well to unity, and scales
     * other defaults appropriately.
     */
    public static void makeLJDefaults() {
        ATOM_SIZE = 1.0;
        ATOM_MASS = 1.0;
        POTENTIAL_WELL = 1.0;
        TEMPERATURE = 1.0;
//        PRESSURE = 1.0;
        TIME_STEP = 0.04;
        BOX_SIZE = 10.0;
        etomica.units.BaseUnit.Length.Sim.TO_PIXELS = 30;
    }
    
    /**
     * Default block size used for error estimation in simulation averages.
     */
    public static int BLOCK_SIZE = 1000;
 
/*    private static String getWorkingDirectory(){
        String dir = System.getProperty("user.dir");
        System.out.println("working directory, in Default: "+dir);
        if(dir.indexOf("VisualCafe") != -1) return "D:\\etomica";
        dir = dir.replace('\\', '/');
        return dir+"/";
    }*/
    
    public static final Parameter.Size SIZE_PARAMETER = new Parameter.Size() {
        double sigma = Default.ATOM_SIZE;
        public double getSigma() {return sigma;}
        public void setSigma(double s) {sigma = s;}
        public Dimension getSigmaDimension() {return Dimension.LENGTH;}
    };
    
    public static final Parameter.Energy ENERGY_PARAMETER = new Parameter.Energy() {
        double epsilon = Default.POTENTIAL_WELL;
        public double getEpsilon() {return epsilon;}
        public void setEpsilon(double e) {epsilon = e;}
        public Dimension getEpsilonDimension() {return Dimension.ENERGY;}
    };
        
    
    public static final Parameter.Mass MASS_PARAMETER = new Parameter.Mass() {
        double mass = Default.ATOM_MASS;
        public double getMass() {return mass;}
        public void setMass(double m) {mass = m;}
        public Dimension getMassDimension() {return Dimension.MASS;}
    };

	//default unit system for I/O (internal calculations are all done in simulation units)
	public static UnitSystem UNIT_SYSTEM = new UnitSystem.Sim();
        
}