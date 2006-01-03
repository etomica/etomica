package etomica.util;

import etomica.units.Bar;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.units.systems.UnitSystem;

/**
 * Class holding fields that define many of the default values used in building
 * a simulation.  
 * 
 * @author kofke
 */

public class Default implements java.io.Serializable {
    
    public double atomSize = 3.0;  //Angstroms
    
    public double atomMass = 40.0; //Daltons
    
    public int moleculeCount = 0;
    
    public double boxSize = 30.0;  //Angstroms
    
    public double temperature = Kelvin.UNIT.toSim(300.);
    
    public double pressure = Bar.UNIT.toSim(1.0);
    
    public double potentialWell = Kelvin.UNIT.toSim(300.);
    
    public double potentialCutoffFactor = 2.5; //dimensionless multiplier for cutoff distance of potential
    
    public double timeStep = 0.05;  //picoseconds 
    
    public int historyPeriod = 100;
    
    public boolean isGraphic = false;
    
    public boolean ignoreOverlap = false;
    
    public Pixel pixelUnit = new Pixel();
    
    /**
     * Default value for doSleep field in ActivityIntegrate class. The default defined
     * here is <tt>false</tt>, indicating that Integrator should not pause during
     * integration loop.  For simulations that use interactive graphics it may
     * be helpful to make doSleep <tt>true</tt> to ensure a responsive
     * interface.  Change in this default has no effect on any classes
     * instantiated prior to the change.
     */
    public boolean doSleep = false;
 
    /**
     * Default block size used for error estimation in simulation averages.
     */
    public int blockSize = 1000;
 
    //default unit system for I/O (internal calculations are all done in simulation units)
    public static UnitSystem UNIT_SYSTEM = new UnitSystem.Sim();
        
    /**
     * Integer array indicating the maximum number of atoms at each depth in the
     * atom hierarchy.  Maximum depth is given by the size of the array.  Each
     * element of array is log2 of the maximum number of child atoms permitted
     * under one atom.  This is used to assign index values to each atom when it
     * is made.  The indexes permit quick comparison of the relative ordering
     * and/or hierarchical relation of any two atoms.  Sum of digits in array
     * should equal 32. For example, {1, 4, 4, 14, 6, 3} indicates 1 SpeciesRoot
     * (this should always be the case), 2^4 = 16 SpeciesMasters (Phases), 2^4 = 16
     * SpeciesAgents (Species), 2^14 =16,384 molecules of a given species in each phase,
     * 2^6 = 64 subgroups per molecule, and 2^3 = 8 atoms per subgroup (note that
     * a molecule need not be structured this way; for example it is ok with this
     * bitLength array to have a monatomic species, there may be 16,384 such monatomic
     * molecules in one phase).
     * 
     * The default is {1, 4, 4, 16, 6, 3}.
     */
    // powers of 2, for reference:
    //  n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |  14  |  15  |  16  |  17   |  18   |  19   |   20    |
    // 2^n| 2 | 4 | 8 | 16| 32| 64|128|256|512|1024|2048|4096|8192|16,384|32,768|65,536|131,072|262,144|524,288|1,048,576|
    // {speciesRoot, phases, species, molecules, groups, atoms}
    public static int[] BIT_LENGTH = new int[] {1, 4, 4, 14, 6, 3};

    
    /**
     * Sets defaults as follows:
     * <ul>
     * <li>atomSize = 1.0
     * <li>atomMass = 1.0
     * <li>potentialWell = 1.0
     * <li>temperature = 1.0
     * <li>pressure = 1.0
     * <li>timeStep = 0.04
     * <li>boxSize = 10
     * <li>pixelUnit = new Pixel(30);
     * </ul>
     */
    public void makeLJDefaults() {
        atomSize = 1.0;
        atomMass = 1.0;
        potentialWell = 1.0;
        temperature = 1.0;
        pressure = 1.0;
        timeStep = 0.04;
        boxSize = 10.0;
        pixelUnit = new Pixel(30);
    }
    
}
