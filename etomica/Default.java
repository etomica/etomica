package etomica;

import etomica.units.*;

public final class Default {
    
    public static String WORKING_DIRECTORY = getWorkingDirectory();
  //  public static String WORKING_DIRECTORY = "D:\\Etomica\\";
//   public static String WORKING_DIRECTORY = "";
    
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
    
    /**
     * Default block size used for error estimation in simulation averages.
     */
    public static int BLOCK_SIZE = 1000;
 
    public static String getWorkingDirectory(){
        String dir = System.getProperty("user.dir");
        System.out.println("working directory, in Default: "+dir);
        if(dir.indexOf("VisualCafe") != -1) return "D:\\etomica";
        dir = dir.replace('\\', '/');
        return dir+"/";
    }
}