package etomica;
import etomica.units.*;

/**
 * Permits selection of temperature from a discrete set of values.  Also has option
 * to select adiabatic.  Enforces selection by manipulating properties of integrator.
 * By default, last integrator added to simulation is the only one controlled by this device.
 *
 */
public class DeviceThermoSelector extends Device implements EtomicaElement {
    
    public String getVersion() {return "DeviceThermoSelector:01.03.23/"+Device.VERSION;}
    /**
     * Descriptive text label to be displayed with the value
     */
    private javax.swing.JComboBox selector;
    private final String adiabaticString = "Adiabatic";
    private Integrator integrator;
    
    private boolean includeAdiabatic = true;
    
    public DeviceThermoSelector() {
        this(Simulation.instance);
    }
    public DeviceThermoSelector(Simulation sim) {
        super(sim);
        selector = new javax.swing.JComboBox();
        setTemperatures(new double[] {200.0, 400.0, 600.0});
        selector.setEditable(false);
        setUnit(new Unit(Kelvin.UNIT));
        
        
/*        selector.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent e){
                Object item = ((javax.swing.JComboBox)e.getSource()).getSelectedItem();
		        if(item==adiabaticString) integrator.setIsothermal(false);
		        else {
		            integrator.setIsothermal(true);
		            integrator.setTemperature(unit.toSim(((Double)item).doubleValue()));
		        }
            }
        });*/
        selector.addItemListener(new java.awt.event.ItemListener() {
		    public void itemStateChanged(java.awt.event.ItemEvent event) {
		        if(integrator==null) return;
		        Object item = event.getItem();
		        if(item==adiabaticString) integrator.setIsothermal(false);
		        else {
		            integrator.setIsothermal(true);
		            integrator.setTemperature(unit.toSim(((Double)item).doubleValue()));
		        }
		    }
		});//end of addItemListener
		
		//add mediator so that last integrator added to simulation is controlled by this device
        sim.elementCoordinator.addMediatorPair(new Mediator.DeviceIntegrator(sim.elementCoordinator) {
            public void add(Integrator integrator) {
                DeviceThermoSelector.this.integrator = integrator;
            }
            public void add(Device device) {
                if(device != DeviceThermoSelector.this) return;
                for(java.util.Iterator ip=mediator.parentSimulation().integratorList.iterator(); ip.hasNext(); ) {
                    Integrator integrator = (Integrator)ip.next();
                    if(integrator.wasAdded())  {//will make last integrator the one
                        DeviceThermoSelector.this.integrator = integrator;
                    }
                }
            }
        });
		
    }//end of constructor
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Select isothermal at a set of temperatures, or adiabatic");
        return info;
    }
    
    public javax.swing.JComboBox getSelector() {return selector;}
    
    public final void setIncludeAdiabatic(boolean b) {
        if(selector.getItemCount() == 0) includeAdiabatic = false;
        if(includeAdiabatic && !b) selector.removeItem(adiabaticString);//remove if was there and removing
        else if(!includeAdiabatic && b) selector.addItem(adiabaticString);//add if wasn't there and adding
        includeAdiabatic = b;
    }
    public boolean isIncludeAdiabatic() {return includeAdiabatic;}
    
    /**
     * Sets the i-th item (counting from 0) in the list as the one currently selected.
     * Argument outside bounds is interpreted as first or last item on list.
     */
    public void setSelected(int i) {
        if(i < 0) i = 0;
        if(i >= selector.getItemCount()) i = selector.getItemCount()-1;
        selector.setSelectedIndex(i);
    }
    
    public void setTemperatures(double[] t) {
        selector.removeAllItems();
        setIncludeAdiabatic(includeAdiabatic);
        for(int i=0; i<t.length; i++) {
            selector.addItem(new Double(t[i]));
        }
        selector.setEnabled(true);
    }

    /**
     * Override superclass setUnit method to update label when unit is changed
     */
    public void setUnit(Unit u) {
        super.setUnit(u);
    }
    
    
    /**
     * Returns the GUI element for display in the simulation.
     * Consists of a combo box used for the selector.
     */
    public java.awt.Component graphic(Object obj) {
        javax.swing.JPanel panel = new javax.swing.JPanel(new java.awt.BorderLayout(0,1));
        panel.add(new javax.swing.JLabel("Temperature ("+unit.symbol()+")"),java.awt.BorderLayout.NORTH);
        panel.add(selector, java.awt.BorderLayout.SOUTH);
        panel.setBorder(new javax.swing.border.EmptyBorder(3,3,3,3));
        return panel;
       // return selector;
    }
}