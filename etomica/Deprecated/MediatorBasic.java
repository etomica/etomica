package simulate;
import javax.swing.JPanel;
        
/**
* Most elementary element coordinator.
* Only actions performed are to place a species agent from each species in each phase,
* and to construct the array of potentials from the registered Potential1 and Potential2 classes.
* Other ElementCoordinator classes can subclass this one, and call the super.go() method from theirs
* to complete these basic actions.
*/
public class MediatorBasic extends Mediator {
    
    public MediatorBasic() {this(Simulation.instance);}
    public MediatorBasic(Simulation sim) {super(sim);}  //constructor
    
    public void go() {
        if(completed) return;
        completed = true;
                
        //Set up arrays of potentials
        int speciesCount = this.parentSimulation().speciesCount();
        if(parentSimulation().potential1 == null) parentSimulation().potential1 = new Potential1[speciesCount];
        if(parentSimulation().potential2 == null) parentSimulation().potential2 = new Potential2[speciesCount][speciesCount];

        Potential1[] potential1 = parentSimulation().potential1;
        Potential2[][] potential2 = parentSimulation().potential2;
        
        Potential1 p1Null = new P1Null();
        Potential2 p2IdealGas = new P2IdealGas();
        for(int i=0; i<parentSimulation().speciesCount(); i++) {
            if (potential1[i] == null) potential1[i] = p1Null;
            for(int j=0; j<speciesCount; j++) {       
                if(potential2[i][j] == null) potential2[i][j] = p2IdealGas;
            }
        }
        if(!Simulation.inEtomica) {
            for(java.util.Iterator iter=Simulation.potential1List.iterator(); iter.hasNext(); ) {
                Potential1 p1 = (Potential1)iter.next();
                potential1[p1.getSpeciesIndex()] = p1;
            }
            for(java.util.Iterator iter=Simulation.potential2List.iterator(); iter.hasNext(); ) {
                Potential2 p2 = (Potential2)iter.next();
                potential2[p2.getSpecies1Index()][p2.getSpecies2Index()] = p2;
                potential2[p2.getSpecies2Index()][p2.getSpecies1Index()] = p2;
            }
        }
                
        //Add a species agent from each species to each phase
        if(!java.beans.Beans.isDesignTime()) {
        for(java.util.Iterator ip=Simulation.phaseList.iterator(); ip.hasNext(); ) {
            Phase p = (Phase)ip.next();
            for(java.util.Iterator is=Simulation.speciesList.iterator(); is.hasNext(); ) {
                Species species = (Species)is.next();
                p.addSpecies(species.makeAgent(p));
            }
        }
        }
                
        //Process graphical elements
        processGraphicalElements();
    }//end of go method
            
    /**
        * Performs default action in processing graphical elements.  
        * Simply adds the graphic from each element to the Simulation.instance.
        * This method can be overridden in subclasses, but should not be invoked there, since
        * it is called in the superclass go() method.
        */
    public void processGraphicalElements() {
	    final javax.swing.JTabbedPane displayPanel = new javax.swing.JTabbedPane();
	    JPanel displayBoxPanel = new JPanel(new java.awt.GridLayout(0,1));
        JPanel devicePanel = new JPanel(new java.awt.GridLayout(0,1));
        for(java.util.Iterator iter=Simulation.graphicalElementList.iterator(); iter.hasNext(); ) {
            Simulation.GraphicalElement element = (Simulation.GraphicalElement)iter.next();
            java.awt.Component component = element.graphic(null);
            if(element instanceof DisplayBox) {
                displayBoxPanel.add(component);
            }
            else if(element instanceof Display) {
                displayPanel.add(((Display)element).getLabel(),component);
            }
            else {
                devicePanel.add(component);
            }
        }
    //       JPanel leftPanel = new JPanel(new GridLayout(2,1));
    //       leftPanel.add(devicePanel);
    //       leftPanel.add(displayBoxPanel);
    //       Simulation.instance.add(leftPanel);
         
        //workaround for JTabbedPane bug in JDK 1.2
        displayPanel.addChangeListener(
            new javax.swing.event.ChangeListener() {
                public void stateChanged(javax.swing.event.ChangeEvent event) {
                    displayPanel.validate();
                }
        });
         
        this.parentSimulation().add(devicePanel);
        this.parentSimulation().add(displayPanel);
        this.parentSimulation().add(displayBoxPanel);
    }//end of processGraphicalElements method
    
}//end of CoordinatorBasic class
        
