package simulate;

import java.awt.Component;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JTable;
import javax.swing.table.*;
import javax.swing.Box;
import javax.swing.JScrollPane;

/**
 * Presents simulation data in a tabular form.  Data is obtained from meter objects.
 * Meters may be identified to the DisplayTable by passing an array of meters (in
 * the setMeter(Meter[]) method, or one at a time via the addMeter(Meter) method.
 * Alternatively, the table may be configured to display all the values of a particular
 * MeterFunction object, identified via the setMeter(MeterFunction) method.
 * Does not permit display of combination of Meters and MeterFunctions.
 */
 
 //setPhase method needs repairing
public class DisplayTable extends Display implements Meter.MultiUser, MeterFunction.User
{
    public JTable table;
    MyTableData dataSource;
    Meter[] meter = new Meter[0];
    MeterFunction meterFunction;
    double[] x;    //used if showing MeterFunction
    double[] y;
    double[] yErr;
    int nMeters = 0;
    private boolean showingPE = false;  //flag to indicate if PotentialEnergy meter should be displayed
    private boolean showingKE = false;  //same for KineticEnergy meter
    boolean showingFunction = false;
    Box panel;
    JButton resetButton;
    /**
     * Shows only averages and error if true, shows only current values if false
     * Default is <code>true</code>
     */
    private boolean showAverages = true;
    
    public DisplayTable() {
        this(Simulation.instance);
    }
    public DisplayTable(Simulation sim) {this(sim,true);}
    public DisplayTable(Simulation sim, boolean showAvgs) {
        super(sim);
        showAverages = showAvgs;
        showingFunction = false;
        setupTable();
    }
    
    private void setupTable() {
        if(panel != null) remove(panel);
        panel = Box.createVerticalBox();
        panel.setSize(100,150);
        dataSource = new MyTableData(showAverages, showingFunction);   //inner class, defined below
        table = new JTable(dataSource);
        panel.add(new JScrollPane(table));
		if(showAverages) {
		    resetButton = new JButton("Reset averages");
		    resetButton.addActionListener(new ActionListener() {
		    public void actionPerformed(java.awt.event.ActionEvent event) {DisplayTable.this.resetAverages();}
		});
            panel.add(resetButton);
        }
        add(panel);
    }
    
    public void setShowAverages(boolean b) {
        showAverages = b;
        setupTable();
    }
    public boolean isShowAverages() {return showAverages;}
    
    public void resetAverages() {
        if(showingFunction) meterFunction.reset();
        else for(int i=0; i<nMeters; i++) {meter[i].reset();}
    }
    public void setResetVisible(boolean b) {resetButton.setVisible(b);}
    public boolean getResetVisible() {return resetButton.isVisible();}
    public void setShowingKE(boolean b) {showingKE = b;}
    public boolean isShowingKE() {return showingKE;}
    public void setShowingPE(boolean b) {showingPE = b;}
    public boolean isShowingPE() {return showingPE;}
    
    public void addMeter(Meter m) {
        if(m instanceof MeterPotentialEnergy && !showingPE) {return;}
        if(m instanceof MeterKineticEnergy && !showingKE) {return;}
        nMeters++;
        Meter[] temp = new Meter[nMeters];
        for(int i=0; i<meter.length; i++) {
            temp[i] = meter[i];
            if(m == meter[i]) return;  //meter is already in table
        }
        temp[nMeters-1] = m;
        meter = temp;
        showingFunction = false;
        meterFunction = null;
        setupTable();
    }
    
    public void setMeters(Meter[] m) {
        meter = m; 
        showingFunction = false;
        meterFunction = null;
    }
    public Meter[] getMeters() {return meter;}
    
    public void setMeterFunction(MeterFunction m) {
        meterFunction = m;
        showingFunction = true;
        meter = null;
        setupTable();
        doUpdate();
    }
    public MeterFunction getMeterFunction() {return meterFunction;}
        
        //might not be working (depends on whether meterList is being kept in phase)
    public void setPhase(Phase p) {
        java.util.LinkedList meterList = p.getMeterList();
        for(java.util.Iterator iter=meterList.iterator(); iter.hasNext(); ) {
            MeterAbstract m = (MeterAbstract)iter.next();
            if(m instanceof Meter) addMeter((Meter)m);
        }
    }
    
    public Component graphic(Object obj) {return panel;}

    //updates x-y values if showing a MeterFunction; otherwise does nothing
    public void doUpdate() {
        if(showingFunction) {
            x = meterFunction.X();
            if(showAverages) {
                y = meterFunction.average();
                yErr = meterFunction.error();
            }
            else {
                y = meterFunction.currentValue();
            }
        }
    }
    public void repaint() {table.repaint();}
    
    private class MyTableData extends AbstractTableModel {
        
        String[] columnNames;
        Class[] columnClasses;
        
        MyTableData(boolean showAverages, boolean showingFunction) {
            if(showingFunction) {  //showing the x-y values of a MeterFunction
                if(showAverages && meterFunction != null) {
                    columnNames = new String[] {meterFunction.getXLabel(), meterFunction.getLabel()+" Average", "y Error"};
                    columnClasses = new Class[] {Double.class, Double.class, Double.class};
                }
                else if(showAverages) {
                    columnNames = new String[] {"x", "y Average", "y Error"};
                    columnClasses = new Class[] {Double.class, Double.class, Double.class};
                }                   
                else if(meterFunction != null) {
                    columnNames = new String[] {"meterFunction.getXLabel()", meterFunction.getLabel()+"Current"};
                    columnClasses = new Class[] {Double.class, Double.class};
                }
                else {
                    columnNames = new String[] {"x", "Current y"};
                    columnClasses = new Class[] {Double.class, Double.class};
                }                   
            }
            else {  //showing a collection of meters
                if(showAverages) {
                    columnNames = new String[] {"Property", "Average", "Error"};
                    columnClasses = new Class[] {String.class, Double.class, Double.class};
                }
                else {
                    columnNames = new String[] {"Property", "Current"};
                    columnClasses = new Class[] {String.class, Double.class};
                }
            }
        }
        
        public Object getValueAt(int row, int column) {
            if(showingFunction) {
                switch(column) {
                    case 0: return new Double(x[row]);    //x-y values are updated in update() method
                    case 1: return new Double(y[row]);    //because we don't want to call currentValue
                    case 2: return new Double(yErr[row]); //or average for each function entry
                    default: return null;
                }
            }
            else {
                Meter m = meter[row];
                switch(column) {
                    case 0: return m.getLabel();
                    case 1: return new Double(showAverages ? m.average() : m.currentValue());
                    case 2: return new Double(m.error());
                    default: return null;
                }
            }
        }
        
        public int getRowCount() {return showingFunction ? meterFunction.getNPoints() : meter.length;}
        public int getColumnCount() {return showAverages ? 3 : 2;}
                
        public String getColumnName(int column) {return columnNames[column];}
        public Class getColumnClass(int column) {return columnClasses[column];}
    }
    
    /**
     * Demonstrates how this class is implemented.  Shows two tables, one 
     * displaying two meters, and one displaying a meterFunction
     */
    public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        Simulation.makeSimpleSimulation();  //for more general simulations, replace this call with
                                            //construction of the desired pieces of the simulation
        //part that is unique to this demonstration
        Meter meter1 = new MeterPressureHard();
        Meter meter2 = new MeterTemperature();
        MeterFunction rdf = new MeterRDF();
        Phase phase = Simulation.instance.phase(0);
        DisplayTable table = new DisplayTable();
        DisplayTable rdfTable = new DisplayTable();
        //end of unique part
                                            
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
        meter1.setPhase(phase);//must come after go() because phase needs to have integrator for this call
        meter2.setPhase(phase);//must come after go() because phase needs to have integrator for this call
        rdf.setPhase(phase);
//        table.setPhase(phase);
        table.addMeter(meter1);
        table.addMeter(meter2);
        rdfTable.setMeterFunction(rdf);
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }

}