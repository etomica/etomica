package etomica;
import java.beans.*;
//Java2 imports
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.LinkedList;
import etomica.utility.Iterator;

/**
 * A property editor for setting the Boundary object of a Phase.
 */

public class BoundaryEditor extends PropertyEditorSupport {
    
    public static String getVersion() {return "01.01.17.0";}

    LinkedList boundaryMakers = new LinkedList();
    
    public BoundaryEditor() {super();}
    public BoundaryEditor(Object obj) {super(obj);}
    
    public String getAsText() {
        return (getValue() != null) ? ((Space.Boundary)getValue()).type().toString() : "";
    }
    
    public void setAsText(String s) {
        String[] options = boundaries();
        for(Iterator iter=boundaryMakers.iterator(); iter.hasNext();) {
            Space.Boundary.Maker maker = (Space.Boundary.Maker)iter.next();
            Space.Boundary.Type[] types = maker.boundaryTypes();
            for(int i=0; i<types.length; i++) {
                if(s.equals(types[i].toString())) {
                    setValue(maker.makeBoundary(types[i]));
                    return;
                }
            }
        }
    }
    
    public String[] getTags() {return boundaries();}
    
    public String[] boundaries() {
        //Get boundary options from space object
        String[] tags =  etomica.utility.StringUtility.toStringArray(Simulation.instance.space().boundaryTypes());
        boundaryMakers.clear();
        boundaryMakers.add(Simulation.instance.space());
        //Loop through all elements in simulation and get any boundary options they provide
        LinkedList elements = Simulation.instance.allElements();
        for(Iterator iter=elements.iterator(); iter.hasNext();) {
            Object obj = iter.next();
            if(obj instanceof Space.Boundary.Maker) {
                Space.Boundary.Maker maker = (Space.Boundary.Maker)obj;
                boundaryMakers.add(maker);
                Space.Boundary.Type[] types = maker.boundaryTypes();
                tags = etomica.utility.StringUtility.arrayCollect(tags, etomica.utility.StringUtility.toStringArray(types));
            }
        }
        return tags;
    }//end of boundaries
    
}