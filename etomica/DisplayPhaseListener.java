package simulate;
import java.util.EventListener;
import java.awt.event.MouseEvent;
import java.awt.event.InputEvent;

/**
 * Interface for a class that listens for events coming from a DisplayPhase.
 * Events of this type usually involve the selection of an atom, molecule, or point in 
 * the phase using the mouse.
 *
 * @see DisplayPhase
 * @see DisplayPhaseEvent
 * @see DisplayPhaseListener.AtomActionWrapper
 */

public interface DisplayPhaseListener extends EventListener {
    
    public void displayPhaseAction(DisplayPhaseEvent dpe);
    
    //end of definition of interface
    
    /** 
     * A wrapper class that takes atom-selection DisplayPhase events and passes them on to an AtomAction
     * An atom-selection MOUSE_PRESSED event (as determined by the DisplayPhase) causes the actionPerformed
     * method of that AtomAction to be called.  This wrapper may be configured to have a MOUSE_RELEASED
     * event or a right-click invoke the retractAction method of the AtomAction (if it implements the
     * <code>Retractable</code> interface).
     * @see AtomAction
     */
     //A main method demonstrating the use of this wrapper is included in AtomAction
    public static class AtomActionWrapper implements DisplayPhaseListener {
        private AtomAction atomAction;
        private boolean retractOnRelease = false;
        private boolean retractOnRightClick = false;
        public AtomActionWrapper(AtomAction aa) {atomAction = aa;}
        /**
         * Accessor method for flag that causes retractAction method to be called on release of the mouse button
         */
        public void setRetractOnRelease(boolean b) {
            if(atomAction instanceof Action.Retractable) retractOnRelease = b;
        }
        /**
         * Accessor method for flag that causes retractAction method to be called on release of the mouse button
         */
        public boolean isRetractOnRelease() {return retractOnRelease;}
        /**
         * Accessor method for flag that causes retractAction method to be called with a right-button press.
         */
        public void setRetractOnRightClick(boolean b) {
            if(atomAction instanceof Action.Retractable) retractOnRightClick = b;
        }
        /**
         * Accessor method for flag that causes retractAction method to be called with a right-button press.
         */
        public boolean isRetractOnRightClick() {return retractOnRightClick;}
        /**
         * DisplayPhaseListener interface method.
         * Interprets mouse event and calls actionPerformed/retractAction method of AtomAction if appropriate
         */
        public void displayPhaseAction(DisplayPhaseEvent dpe) {
            MouseEvent mouseEvent = dpe.getMouseEvent();
            if(mouseEvent == null) return;
            Atom atom = dpe.getAtom();
            if(atom == null) return;
            switch(mouseEvent.getID()) {
                case MouseEvent.MOUSE_PRESSED:
                    if(retractOnRightClick && ((mouseEvent.getModifiers() & InputEvent.BUTTON3_MASK) != 0)) 
                            ((Action.Retractable)atomAction).retractAction();
                    else atomAction.actionPerformed(atom);
                    break;
                case MouseEvent.MOUSE_RELEASED:
                    if(retractOnRelease) ((Action.Retractable)atomAction).retractAction();
                    break;
                default:
                    break;
            }//end of switch statement
        }//end of displayPhaseAction
    }//end of AtomActionWrapper

}//end of DisplayPhaseListener