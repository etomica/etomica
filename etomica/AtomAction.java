// This class includes a main method demonstrating its use.
package etomica;
import java.awt.Color;  //for the color-change action
//imports for main method
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.Frame;

/**
 * Base class for classes that perform some elementary action on an atom.
 * These classes see use in several ways:
 * <ul>
 * <li>they can be passed to the allAtoms method of atom iterators, which then performs
 * the specified action on all the atoms of the iterator.  
 * <li>they can be used in a DisplayPhase.AtomActionWrapper, to specify some action in response to 
 * selection of an atom by the mouse.  
 * <li>they may be used to generate a Monte Carlo trial in an MCMove object.
 * 
 * @author David Kofke
 * @see Atom.Iterator
 * @see DisplayPhase.AtomActionWrapper
 */

public abstract class AtomAction extends etomica.Action {
    
    public static String getVersion() {return "01.01.17.0/"+Action.getVersion();}

    protected Atom atom;
    public void setAtom(Atom a) {atom = a;}
    public Atom getAtom() {return atom;}

    /**
     * Performs the defined action using the atom most recently specified by setAtom or by the last call to actionPerformed(Atom a).
     * Performs no action if the atom is null.
     */
    public void actionPerformed() {if(atom != null) actionPerformed(atom);}
    
    /**
     * Method that defines the action to be performed on the atom
     * @param a Atom passed to method by iterator
     */
    public abstract void actionPerformed(Atom a);
        
    //***** end of Action methods; begin definition of subclasses *****//

    /**
     * Changes the color of an atom to some specified value.
     */
    public static class ChangeColor extends AtomAction implements Action.Retractable {
        private Color newColor;
        private Color oldColor;
        public ChangeColor() {this(Color.red);}
        public ChangeColor(Color c) {setNewColor(c);}
        public void setNewColor(Color c) {newColor = c;}
        public Color getNewColor() {return newColor;}
        public void actionPerformed(Atom a) {
            atom = a;
            oldColor = a.getColor();
            a.setColor(newColor);
        }
        public void retractAction() {
            if(atom != null) atom.setColor(oldColor);
        }
    }
    
    /**
     * Translates the atom by the amount it would move in free (ballistic) flight for a specified time interval.
     * Uses the atom's current momentum to determine this displacement.
     */
    public static class FreeFlight extends AtomAction {
        private double tStep = 0.0;
        public void actionPerformed(Atom a) {
            if(a.isStationary()) {return;}  //skip if atom is stationary
            a.r.PEa1Tv1(tStep*a.rm(),a.p);  // r += tStep*p/m
        }
        public void actionPerformed(Atom a, double t) {
            tStep = t;
            actionPerformed(a);
        }
        public void setTStep(double t) {tStep = t;}
        public double getTStep() {return tStep;}
    } //end of FreeFlight
        
    public static class Translate extends AtomAction {
        protected Space.Vector displacement;
            
        public Translate(Space space) {
            super();
            displacement = space.makeVector();
        }
            
        public final void actionPerformed(Atom a) {a.r.PE(displacement);}
        public void actionPerformed(Atom a, Space.Vector d) {a.r.PE(d);}
        public final void setDisplacement(Space.Vector d) {displacement.E(d);}
    }//end of Translate
    
//        public static class Displace extends Translate {
//            Atom atom;
//            public void actionPerformed(Atom a) {atom = a; a.r.PE(displacement);}
//            public void actionPerformed(Atom a, Space.Vector d) {
//                displacement.E(d);
//                a.r.PE(displacement);
//            }
//            public void retractAction() {a.r.ME(displacement);}
//        }//end of Displace

    /**
     * Demonstrates how this class can be implemented with a DisplayPhase event using a AtomActionWrapper.
     * Hold the 'a' key while pressing mouse button near an atom to change its color; relase of the 
     * mouse button reverts to original color.
     */
    public static void main(String[] args) {
        final Frame f = new Frame();   //create a window
        f.setSize(600,350);
        
      //make a simple simulation for this example
        Simulation.makeSimpleSimulation();
      //get a handle to the display in the simple simulation
        final DisplayPhase display = (DisplayPhase)Simulation.instance.displayList.get(0); 
      //create an instance of the AtomAction that is being demonstrated here
        AtomAction.ChangeColor colorChanger = new AtomAction.ChangeColor();
      //wrap the action in a DisplayPhaseListener so it responds to MousePressed events
        DisplayPhaseListener.AtomActionWrapper wrapper = new DisplayPhaseListener.AtomActionWrapper(colorChanger);
      //set the wrapper so that release of the mouse button retracts the action
        wrapper.setRetractOnRelease(true);
      //make the wrapper a listener to the DisplayPhase, and add listener to repaint display
        display.addDisplayPhaseListener(wrapper);
        display.addDisplayPhaseListener(new DisplayPhaseListener() {
            public void displayPhaseAction(DisplayPhaseEvent e) {display.repaint();}
        });
      //add the simulation graphic elements to the frame
        f.add(Simulation.instance);
      //display the frame
        f.pack();
        f.show();
        f.addWindowListener(new WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(WindowEvent e) {System.exit(0);}
        });
    }
} //end of AtomAction   