package etomica;

import java.awt.Color;
import etomica.units.*;

public final class Default {
    
    public static String WORKING_DIRECTORY = getWorkingDirectory();
    
    public static String CLASS_DIRECTORY = WORKING_DIRECTORY + "etomica";
    
    public static String IMAGE_DIRECTORY = "file:/" + WORKING_DIRECTORY + "images/";
    
    public static String HELP_FILE = "http://www.ccr.buffalo.edu/etomica/help.html";
        
    public static String JAVADOC_FILE = "http://www.ccr.buffalo.edu/etomica/JavaDoc/index.html";

    public static double ATOM_SIZE = 3.0;  //Angstroms
    
    public static double ATOM_MASS = 40.0; //Daltons
    
    public static Color ATOM_COLOR = Color.black;
    
    public static int MOLECULE_COUNT = 20;
    
    public static double BOX_SIZE = 30.0;  //Angstroms
    
    public static double TEMPERATURE = Kelvin.UNIT.toSim(300.);
    
    public static double PRESSURE = 1.0;
    
    public static double POTENTIAL_WELL = Kelvin.UNIT.toSim(300.);
    
    public static double POTENTIAL_CUTOFF = 2.5; //dimensionless multiplier for cutoff distance of potential
    
    public static double TIME_STEP = 0.05;  //picoseconds
 
    public static String getWorkingDirectory(){
        String dir = System.getProperty("user.dir");
        dir = dir.replace('\\', '/');
        return dir+"/";
    }
}