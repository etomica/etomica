/**
 * EventActions
 *
 * This class provides static action listeners for most major events that can occur during a simulation.
 * These include mouse, keyboard, component, and property changes.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */

package simulate.gui;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class EventActions {
    /**
     * static action listener for property changes
     */
    public static final ActionListener PROPCHANGE = new PropChangeAction();
    
    /**
     * static action listener for dragging of the mouse
     */
    public static final ActionListener MOUSEDRAGGED = new MouseDraggedAction();
    
    /**
     * static action listener for movement of the mouse
     */
    public static final ActionListener MOUSEMOVED = new MouseMovedAction();
    
    /**
     * static action listener for a component gaining focus
     */
    public static final ActionListener FOCUSGAINED = new FocusGainedAction();
    
    /**
     * static action listener for a component losing focus
     */
    public static final ActionListener FOCUSLOST = new FocusLostAction();
    
    /**
     * static action listener for clicking of the mouse
     */
    public static final ActionListener MOUSECLICKED = new MouseClickedAction();
    
    /**
     * static action listener for the entering of the cursor into a component's bounds
     */
    public static final ActionListener MOUSEENTERED = new MouseEnteredAction();
    
    /**
     * static action listener for the exiting of the cursor from a component's bounds
     */
    public static final ActionListener MOUSEEXITED = new MouseExitedAction();
    
    /**
     * static action listener for pressing a mouse button
     */
    public static final ActionListener MOUSEPRESSED = new MousePressedAction();
    
    /**
     * static action listener for releasing a mouse button
     */
    public static final ActionListener MOUSERELEASED = new MouseReleasedAction();
    
    /**
     * static action listener for changing the caret position
     */
    public static final ActionListener CARETPOSITIONCHANGED = new CaretPositionChangedAction();
    
    /**
     * static action listener for changing of the text of a textfield
     */
    public static final ActionListener INPUTMETHODTEXTCHANGED = new InputMethodTextChangedAction();
    
    /**
     * static action listener for adding a component
     */
    public static final ActionListener COMPONENTADDED = new ComponentAddedAction();
    
    /**
     * static action listener for removing a component
     */
    public static final ActionListener COMPONENTREMOVED = new ComponentRemovedAction();
    
    /**
     * static action listener for pressing a key on the keyboard
     */
    public static final ActionListener KEYPRESSED = new KeyPressedAction();
    
    /**
     * static action listener for releasing a key on the keyboard
     */
    public static final ActionListener KEYRELEASED = new KeyReleasedAction();
    
    /**
     * static action listener for typing a key on the keyboard
     */
    public static final ActionListener KEYTYPED = new KeyTypedAction();
    
    /**
     * static action listener for hiding a component
     */
    public static final ActionListener COMPONENTHIDDEN = new ComponentHiddenAction();
    
    /**
     * static action listener for moving a component
     */
    public static final ActionListener COMPONENTMOVED = new ComponentMovedAction();
    
    /**
     * static action listener for resizing a component
     */
    public static final ActionListener COMPONENTRESIZED = new ComponentResizedAction();
    
    /**
     * static action listener for showing a component
     */
    public static final ActionListener COMPONENTSHOWN = new ComponentShownAction();
    
    /**
     * static class that handles the property change event
     */
    private static class PropChangeAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of PropChangeAction

    /**
     * static class that handles the mouse dragged event
     */
    private static class MouseDraggedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of MouseDraggedAction
    
    /**
     * static class that handles the mouse moved event
     */
    private static class MouseMovedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of MouseMovedAction

    /**
     * static class that handles the focus gained event
     */
    private static class FocusGainedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of FocusGainedAction
    
    /**
     * static class that handles the focus lost event
     */
    private static class FocusLostAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of FocusLostAction

    /**
     * static class that handles the mouse clicked event
     */
    private static class MouseClickedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of MouseClickedAction
    
    /**
     * static class that handles the mouse entered event
     */
    private static class MouseEnteredAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of MouseEnteredAction

    /**
     * static class that handles the mouse exited event
     */
    private static class MouseExitedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of MouseExitedAction
    
    /**
     * static class that handles the mouse pressed event
     */
    private static class MousePressedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of MousePressedAction

    /**
     * static class that handles the mouse released event
     */
    private static class MouseReleasedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of MouseReleasedAction
    
    /**
     * static class that handles the caret position changed event
     */
    private static class CaretPositionChangedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of CaretPositionChangedAction

    /**
     * static class that handles the input method text changed event
     */
    private static class InputMethodTextChangedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of InputMethodTextChangedAction
    
    /**
     * static class that handles the component added event
     */
    private static class ComponentAddedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of ComponentAddedAction

    /**
     * static class that handles the component removed event
     */
    private static class ComponentRemovedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of ComponentRemovedAction
    
    /**
     * static class that handles the key pressed event
     */
    private static class KeyPressedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of KeyPressedAction

    /**
     * static class that handles the key released event
     */
    private static class KeyReleasedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of KeyReleasedAction
    
    /**
     * static class that handles the key typed event
     */
    private static class KeyTypedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of KeyTypedAction

    /**
     * static class that handles the component hidden event
     */
    private static class ComponentHiddenAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of ComponentHiddenAction
    
    /**
     * static class that handles the component moved event
     */
    private static class ComponentMovedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of ComponentMovedAction

    /**
     * static class that handles the component resized event
     */
    private static class ComponentResizedAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of ComponentResizedAction
    
    /**
     * static class that handles the component shown event
     */
    private static class ComponentShownAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of ComponentShownAction
}// end of EventActions class
    