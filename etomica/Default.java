package etomica;

import java.awt.Color;
import etomica.units.*;

public final class Default {
    
    public static String WORKING_DIRECTORY = getWorkingDirectory();
 //   public static String WORKING_DIRECTORY = "D:\\Etomica\\";
 //   public static String WORKING_DIRECTORY = "";
    
    public static String CLASS_DIRECTORY = WORKING_DIRECTORY + "etomica";
    
    public static String IMAGE_DIRECTORY = "file:/" + WORKING_DIRECTORY + "images/";
    
    public static String HELP_FILE = "http://www.ccr.buffalo.edu/etomica/help.html";
        
    public static String JAVADOC_FILE = "http://www.ccr.buffalo.edu/etomica/JavaDoc/index.html";

    public static boolean DISPLAY_USE_OPENGL = true;
    
    public static double ATOM_SIZE = 3.0;  //Angstroms
    
    public static double ATOM_MASS = 40.0; //Daltons
    
    public static Color ATOM_COLOR = Color.black;
    
    public static int MOLECULE_COUNT = 20;
    
    public static double BOX_SIZE = 30.0;  //Angstroms
    
    public static double TEMPERATURE = Kelvin.UNIT.toSim(300.);
    
    public static double PRESSURE = Bar.UNIT.toSim(1.0);
    
    public static double POTENTIAL_WELL = Kelvin.UNIT.toSim(300.);
    
    public static double POTENTIAL_CUTOFF_FACTOR = 2.5; //dimensionless multiplier for cutoff distance of potential
    
    public static double TIME_STEP = 0.05;  //picoseconds 
    
    public static int HISTORY_PERIOD = 100;
    
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