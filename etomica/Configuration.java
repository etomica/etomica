package simulate;
import java.util.*;
import java.awt.*;

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
