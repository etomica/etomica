package etomica.gui;
import etomica.Space;

public class SpaceSelectionFrame extends javax.swing.JInternalFrame implements java.awt.event.ActionListener {

    Listener listener; //inner interface defined below
    javax.swing.ButtonGroup dimension;
    
    public SpaceSelectionFrame(Listener s) {
        super("Space",
	                true, // resizable
	                true, // closable
	                true, // maximizable
	                true); // iconifiable
        listener = s;	        
        getContentPane().setLayout(new java.awt.GridLayout(0,1));
	    dimension = new javax.swing.ButtonGroup();
    		    
        //This section creates the radio buttons for all the classes that subclass space.class, makes
        //them mutually exclusive, and adds buttonlisteners so that the selected button is known
	    for(int i=0; i<Etomica.spaceClasses.length; i++) {
            String name = Etomica.spaceClasses[i].getName();
            int idx = 12;//strip off etomica.Space prefix
            name = name.substring(idx+1);
            MyButton button = new MyButton(Etomica.spaceClasses[i],name,i==0);
            dimension.add(button);
            getContentPane().add(button);
        }// end of radio button creation
                
        /**
        * This section creates the "OK" button that when pressed creates an instance of 
        * simulation.instance, adds it to an internal frame, and opens the corresponding editor 
        * window. 
        */
        javax.swing.JButton oK = new javax.swing.JButton("OK");
        oK.addActionListener(this);
        
        getContentPane().add(oK);
        reshape(412, 200, 200, 200);
        setVisible(true);

		Etomica.DesktopFrame.desktop.add(this);
        try{setSelected(true);}
        catch(java.beans.PropertyVetoException e){} // attempt was vetoed

    }//end of constructor
    
    /**
     * Simple extension of JRadioButton to hold its corresponding space class object.
     */
    private static class MyButton extends javax.swing.JRadioButton {
        Class spaceClass;
        MyButton(Class c, String name, boolean selected) {
            super(name, selected);
            spaceClass = c;
        }
    }
            
    /**
     * Activities performed when the space selection is completed by pressing the "OK" button.
     */
    public void actionPerformed(java.awt.event.ActionEvent e){
        Space space = null;
        MyButton button = null;
        try{
  //          MyButton button = (MyButton)(dimension.getSelection().getSelectedObjects()[0]);
            java.awt.Component[] allComponents = getContentPane().getComponents();
            // Run through all the space radio buttons
            for(int j = 0; j < allComponents.length-1; j++){
                // See which one is selected
                if (((javax.swing.JRadioButton)allComponents[j]).isSelected()){
                // Create SimulationFrame
                    button = (MyButton)allComponents[j];
                }// end of manipulations necessary based on the selected space class
            }// end of loop that runs through all of the space classes
            
            space = (Space)button.spaceClass.newInstance();
        }
        catch(IllegalAccessException ex){}
        catch(InstantiationException ex) {}
        try {this.setClosed(true);}
        catch(java.beans.PropertyVetoException ex) {}
        
        listener.spaceSelectionAction(space);
    }// end of actionPerformed
    
    /**
     * Interface describing the class that created this space selection frame.
     * It has a method that is invoked when the selection frame selects and
     * instantiates a new space.
     */
    public interface Listener {
        public void spaceSelectionAction(Space s);
    }
}//end of SpaceSelectionFrame
        