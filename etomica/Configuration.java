package simulate;
import java.util.*;
import java.awt.*;

/**
 * General class for assignment of coordinates to all molecules/atoms in phase
 */

public abstract class Configuration extends Component{

    Vector species = new Vector();
    
    public Configuration(){
    }
    
    public Configuration(Species s){
        species.addElement(s);
        initializeCoordinates();
    }
    
    public void add(Species s){
        if(s instanceof SpeciesWalls) {return;}
        species.addElement(s);
        initializeCoordinates();
    }

    public abstract void initializeCoordinates();
    	
}
