/**
 * EtomicaMenuBar
 *
 * This class creates the main menu bar for the Etomica simulation environment.  This includes, giving
 * titles, adding JMenu and JMenuItem instances, adding hotkeys, and adding actionListeners that link to
 * static action listeners defined in the <JMenu name> + Actions classes.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */

package simulate.gui;

import simulate.Simulation;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComponent;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JMenu;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
import javax.swing.plaf.ComponentUI;
import javax.swing.plaf.basic.BasicToolBarUI;

public class EtomicaMenuBar extends JMenuBar {
 
	/**
	 * Menu Handles
	 */
	static final JMenuBar JMenuBar1 = new JMenuBar();
	static final JMenu fileMenu = new JMenu("File");
	static final JMenu editMenu = new JMenu("Edit");
	static final JMenu eventMenu = new JMenu("Event");
	static final JMenu viewMenu = new JMenu("View");
	static final JMenu windowMenu = new JMenu("Window");
	static final JMenu simulationMenu = new JMenu("Simulation");
	static final JMenu helpMenu = new JMenu("Help");
	
	/**
	 * File MenuItem Handles
	 */
	static final JMenuItem newSimulationItem = new JMenuItem("New Simulation");
	static final JMenuItem serAppletItem = new JMenuItem("Serialize (Applet Form)");
	static final JMenuItem serEditItem = new JMenuItem("Serialize (Edit Form)");
	static final JMenuItem loadItem = new JMenuItem("Load");
	static final JMenuItem printItem = new JMenuItem("Print");
	static final JMenuItem clearItem = new JMenuItem("Clear");
	static final JMenuItem exitItem = new JMenuItem("Exit");
	static final JSeparator JSeparator1 = new JSeparator();
	
	/**
	 * Edit MenuItem Handles    
	 */
	static final JMenuItem cutItem = new JMenuItem("Cut");
	static final JMenuItem copyItem = new JMenuItem("Copy");
	static final JMenuItem pasteItem = new JMenuItem("Paste");
	static final JMenuItem customizeItem = new JMenuItem("Customize");
	static final JMenuItem reportItem = new JMenuItem("Report");
	static final JMenuItem bindPropertyItem = new JMenuItem("Bind Property");
    
    /**
     * View MenuItem Handles
	 */
	static final JMenuItem workBookItem = new JMenuItem("Workbook");
	static final JMenuItem propertyListItem = new JMenuItem("Property List");
	static final JMenuItem componentLibraryItem = new JMenuItem("Component Library");
	static final JMenuItem classBrowserItem = new JMenuItem("Class Browser");
	static final JMenuItem messagesItem = new JMenuItem("Messages");
   
    /**
     * Event MenuItem Handles
     */
    static final JMenu propChangeMenu = new JMenu("Property Change");
    static final JMenu mouseMotionMenu = new JMenu("Mouse Motion");
    static final JMenu focusMenu = new JMenu("Focus");
    static final JMenu mouseMenu = new JMenu("Mouse");
    static final JMenu inputMethodMenu = new JMenu("Input Method");
    static final JMenu containerMenu = new JMenu("Container");
    static final JMenu keyMenu = new JMenu("Key");
    static final JMenu componentMenu = new JMenu("Component");
    
    /**
     * Event SubMenuItem Handles
     */
    static final JMenuItem propChangeItem = new JMenuItem("Property Changed");
    static final JMenuItem mouseDraggedItem = new JMenuItem("Mouse Dragged");
    static final JMenuItem mouseMovedItem = new JMenuItem("Mouse Moved");
    static final JMenuItem focusGainedItem = new JMenuItem("Focus Gained");
    static final JMenuItem focusLostItem = new JMenuItem("Focus Lost");
    static final JMenuItem mouseClickedItem = new JMenuItem("Mouse Clicked");
    static final JMenuItem mouseEnteredItem = new JMenuItem("Mouse Entered");
    static final JMenuItem mouseExitedItem = new JMenuItem("Mouse Exited");
    static final JMenuItem mousePressedItem = new JMenuItem("Mouse Pressed");
    static final JMenuItem mouseReleasedItem = new JMenuItem("Mouse Released");
    static final JMenuItem caretPositionChangedItem = new JMenuItem("Caret Position Changed");
    static final JMenuItem inputMethodTextChangedItem = new JMenuItem("Input Method Text Changed");
    static final JMenuItem componentAddedItem = new JMenuItem("Component Added");
    static final JMenuItem componentRemovedItem = new JMenuItem("Component Removed");
    static final JMenuItem keyPressedItem = new JMenuItem("Key Pressed");
    static final JMenuItem keyReleasedItem = new JMenuItem("Key Released");
    static final JMenuItem keyTypedItem = new JMenuItem("Key Typed");
    static final JMenuItem componentHiddenItem = new JMenuItem("Component Hidden");
    static final JMenuItem componentMovedItem = new JMenuItem("Component Moved");
    static final JMenuItem componentResizedItem = new JMenuItem("Component Resized");
    static final JMenuItem componentShownItem = new JMenuItem("Component Shown");
    
    /**
     * SimulationItem Menu Handles
     */
    static final JMenuItem selectSpaceItem = new JMenuItem("Select Space");
    static final JMenuItem editSimulationItem = new JMenuItem("Edit Simulation");
    
    /**
     * Window MenuItem Handles
     */
    static final JMenuItem nextItem = new JMenuItem("Next");
	static final JMenuItem cascadeItem = new JMenuItem("Cascade");
	static final JMenuItem tileItem = new JMenuItem("Tile");    
    static final JMenuItem dragOutlineItem = new JMenuItem("Drag Outline");   
    
    /**
     * Help MenuItem Handles
     */
    static final JMenuItem helpItem = new JMenuItem("Help");
    static final JMenuItem aboutItem = new JMenuItem("About");
    
    /**
     * Constructor that adds all new JMenus to the JMenuBar, then adds all JMenuItems to each JMenu.
     * ActionListeners are also added to each JMenuItem, as well as hotkeys.
     */
    public EtomicaMenuBar() {
//        BasicToolBarUI btbUI = new BasicToolBarUI();
//       this.setUI((ComponentUI)btbUI);
     
        /**
         * File Menu naming, adding to menu, and listener creating 
         */
        this.add(fileMenu);
        newSimulationItem.addActionListener(FileActions.NEW_SIMULATION);
        fileMenu.add(newSimulationItem);
        serEditItem.addActionListener(FileActions.SEREDIT);
		serEditItem.setAccelerator(KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_S, java.awt.Event.CTRL_MASK));
        serEditItem.setEnabled(false);
        fileMenu.add(serEditItem);
        serAppletItem.addActionListener(FileActions.SERAPPLET);
        serAppletItem.setEnabled(false);
        fileMenu.add(serAppletItem);
		loadItem.setActionCommand("Load");
		loadItem.setAccelerator(KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_L, java.awt.Event.CTRL_MASK));
		loadItem.setMnemonic((int)'L');
        loadItem.addActionListener(FileActions.LOAD);
        fileMenu.add(loadItem);
        printItem.setActionCommand("Print");
		printItem.setAccelerator(KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_P, java.awt.Event.CTRL_MASK));
		printItem.setMnemonic((int)'P');
        printItem.addActionListener(FileActions.PRINT);
        fileMenu.add(printItem);
        clearItem.addActionListener(FileActions.CLEAR);
        fileMenu.add(clearItem);
        exitItem.setActionCommand("Exit");
		exitItem.setAccelerator(KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_F4, java.awt.Event.ALT_MASK));
        exitItem.addActionListener(FileActions.EXIT);
        fileMenu.add(exitItem);
        
        /**
         * Edit Menu naming, adding to menu, and listener creating 
         */
        this.add(editMenu);
        cutItem.setActionCommand("Cut");
		cutItem.setAccelerator(KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_X, java.awt.Event.CTRL_MASK));
        cutItem.addActionListener(EditActions.CUT);
        editMenu.add(cutItem);
        copyItem.setActionCommand("Copy");
		copyItem.setAccelerator(KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_C, java.awt.Event.CTRL_MASK));
		copyItem.setMnemonic((int)'C');
        copyItem.addActionListener(EditActions.COPY);
        editMenu.add(copyItem);
        pasteItem.setActionCommand("Paste");
		pasteItem.setAccelerator(KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_V, java.awt.Event.CTRL_MASK));
        pasteItem.addActionListener(EditActions.PASTE);
        editMenu.add(pasteItem);
        customizeItem.addActionListener(EditActions.CUSTOMIZE);
        editMenu.add(customizeItem);
        reportItem.addActionListener(EditActions.REPORT);
        editMenu.add(reportItem);
        
        /**
         * View Menu naming, adding to menu, and listener creating 
         */
        this.add(viewMenu);
        workBookItem.addActionListener(ViewActions.WORKBOOK);
        workBookItem.setActionCommand("Workbook");
        viewMenu.add(workBookItem);
        propertyListItem.addActionListener(ViewActions.PROPERTYLIST);
        propertyListItem.setActionCommand("PropertyList");
		propertyListItem.setAccelerator(KeyStroke.getKeyStroke("F4"));
		propertyListItem.setMnemonic((int)'P');
        viewMenu.add(propertyListItem);
        componentLibraryItem.addActionListener(ViewActions.COMPONENTLIBRARY);
        componentLibraryItem.setActionCommand("Component Library");
		componentLibraryItem.setMnemonic((int)'L');
        viewMenu.add(componentLibraryItem);
        classBrowserItem.addActionListener(ViewActions.CLASSBROWSER);
        classBrowserItem.setActionCommand("Class Browser");
		classBrowserItem.setAccelerator(KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_C, java.awt.Event.CTRL_MASK + java.awt.Event.SHIFT_MASK));
		classBrowserItem.setMnemonic((int)'C');
        viewMenu.add(classBrowserItem);
        messagesItem.addActionListener(ViewActions.MESSAGES);
        messagesItem.setActionCommand("Messages");
		messagesItem.setAccelerator(KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_M, java.awt.Event.CTRL_MASK + java.awt.Event.SHIFT_MASK));
		messagesItem.setMnemonic((int)'M');
        viewMenu.add(messagesItem);
        
        /**
         *  Event Menu listed under Edit Menu naming, adding to menu, and listener creating 
         */
        editMenu.add(eventMenu);
        
        /**
         * Property Change Menu (Sub Menu of Event Menu) naming, adding to menu, and listener creating 
         */
        eventMenu.add(propChangeMenu);
        propChangeItem.addActionListener(EventActions.PROPCHANGE);
        propChangeMenu.add(propChangeItem);
        
        /**
         * Mouse Motion Menu (Sub Menu of Event Menu) naming, adding to menu, and listener creating 
         */
        eventMenu.add(mouseMotionMenu);
        mouseDraggedItem.addActionListener(EventActions.MOUSEDRAGGED);
        mouseMotionMenu.add(mouseDraggedItem);
        mouseMovedItem.addActionListener(EventActions.MOUSEMOVED);
        mouseMotionMenu.add(mouseMovedItem);

        /**
         * Focus Menu (Sub Menu of Event Menu) naming, adding to menu, and listener creating 
         */
        eventMenu.add(focusMenu);
        focusGainedItem.addActionListener(EventActions.FOCUSGAINED);
        focusMenu.add(focusGainedItem);
        focusLostItem.addActionListener(EventActions.FOCUSLOST);
        focusMenu.add(focusLostItem);

        /**
         * Mouse Menu (Sub Menu of Event Menu) naming, adding to menu, and listener creating 
         */
        eventMenu.add(mouseMenu);
        mouseClickedItem.addActionListener(EventActions.MOUSECLICKED);
        mouseMenu.add(mouseClickedItem);
        mouseEnteredItem.addActionListener(EventActions.MOUSEENTERED);
        mouseMenu.add(mouseEnteredItem);
        mouseExitedItem.addActionListener(EventActions.MOUSEEXITED);
        mouseMenu.add(mouseExitedItem);
        mousePressedItem.addActionListener(EventActions.MOUSEPRESSED);
        mouseMenu.add(mousePressedItem);
        mouseReleasedItem.addActionListener(EventActions.MOUSERELEASED);
        mouseMenu.add(mouseReleasedItem);
 
        /**
         * Input Method Menu (Sub Menu of Event Menu) naming, adding to menu, and listener creating 
         */
        eventMenu.add(inputMethodMenu);
        caretPositionChangedItem.addActionListener(EventActions.CARETPOSITIONCHANGED);
        inputMethodMenu.add(caretPositionChangedItem);
        inputMethodTextChangedItem.addActionListener(EventActions.INPUTMETHODTEXTCHANGED);
        inputMethodMenu.add(inputMethodTextChangedItem);
 
        /**
         * Container Menu (Sub Menu of Event Menu) naming, adding to menu, and listener creating 
         */
        eventMenu.add(containerMenu);
        componentAddedItem.addActionListener(EventActions.COMPONENTADDED);
        containerMenu.add(componentAddedItem);
        componentRemovedItem.addActionListener(EventActions.COMPONENTREMOVED);
        containerMenu.add(componentRemovedItem);

        /**
         * Key Menu (Sub Menu of Event Menu) naming, adding to menu, and listener creating 
         */
        eventMenu.add(keyMenu);
        keyPressedItem.addActionListener(EventActions.KEYPRESSED);
        keyMenu.add(keyPressedItem);
        keyReleasedItem.addActionListener(EventActions.KEYRELEASED);
        keyMenu.add(keyReleasedItem);
        keyTypedItem.addActionListener(EventActions.KEYTYPED);
        keyMenu.add(keyTypedItem);
 
        /**
         * Component Menu (Sub Menu of Event Menu) naming, adding to menu, and listener creating 
         */
        eventMenu.add(componentMenu);
        componentHiddenItem.addActionListener(EventActions.COMPONENTHIDDEN);
        componentMenu.add(componentHiddenItem);
        componentMovedItem.addActionListener(EventActions.COMPONENTMOVED);
        componentMenu.add(componentMovedItem);
        componentResizedItem.addActionListener(EventActions.COMPONENTRESIZED);
        componentMenu.add(componentResizedItem);
        componentShownItem.addActionListener(EventActions.COMPONENTSHOWN);
        componentMenu.add(componentShownItem);
        bindPropertyItem.addActionListener(EditActions.BINDPROPERTY);
        editMenu.add(bindPropertyItem);
                
        /**
         * Window Menu naming, adding to menu, and listener creating 
         */
        this.add(windowMenu);
        nextItem.addActionListener(WindowsActions.NEXT);
        windowMenu.add(nextItem);
        cascadeItem.addActionListener(WindowsActions.CASCADE);
        windowMenu.add(cascadeItem);
        tileItem.addActionListener(WindowsActions.TILE);
        windowMenu.add(tileItem);
//        dragOutlineItem.addActionListener(WindowsActions.DRAG);
//        windowMenu.add(dragOutlineItem);
      
        /**
         * Simulate Menu naming, adding to menu, and listener creating 
         */
        this.add(simulationMenu);
        selectSpaceItem.addActionListener(SimulateActions.SELECTSPACE);
        simulationMenu.add(selectSpaceItem);
        editSimulationItem.setEnabled(false);
        editSimulationItem.addActionListener(SimulateActions.EDITSIMULATION);
        simulationMenu.add(editSimulationItem);
        simulationMenu.add(new javax.swing.JSeparator());//instantiated simulations will be placed after this separator, as they are created

        /**
         * Help Menu naming, adding to menu, and listener creating 
         */
        this.add(helpMenu);
		helpItem.setActionCommand("Help");
		helpItem.setMnemonic((int)'H');
		helpItem.addActionListener(HelpActions.HELP);
		helpMenu.add(helpItem);
        aboutItem.addActionListener(HelpActions.ABOUT);
		helpMenu.add(aboutItem);
	}// end of EtmoicaMenuBar constructor
}// end of EtomicaMenuBar class