package etomica;

import java.awt.Component;
import javax.swing.JTable;
import javax.swing.table.*;
import javax.swing.Box;
import javax.swing.JScrollPane;
import etomica.units.Unit;

/**
 * Presents function data in a tabular form.  Data is obtained from DataSource objects.
 * Sources may be identified to the Display by passing an array of DataSource objects (in
 * the setDataSources method, or one at a time via the addDataSource method).
 */
public class DisplayTableFunction extends DisplayDataSources implements EtomicaElement
{
    public String getVersion() {return "DisplayTableFunction:01.05.29/"+super.getVersion();}

    public JTable table;
    Box panel;
    
    private boolean showXColumn = true;
        
    public DisplayTableFunction() {
        this(Simulation.instance);
    }
    public DisplayTableFunction(Simulation sim)  {
        super(sim);
        setupDisplay();
        setLabel("Function");
        setWhichValue(MeterAbstract.AVERAGE);
        table = new JTable();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Tabular display of function data from several sources");
        return info;
    }
    
    /**
     * Implementation of parent's abstract method to set up table whenever 
     * data sources are added or changed.
     */
    protected void setupDisplay() {
        if(panel == null) panel = Box.createVerticalBox();
        else panel.removeAll();
        panel.setSize(100,150);
        if(ySource == null || ySource.length == 0) return;
        MyTableData tableSource = new MyTableData();   //inner class, defined below
        table = new JTable(tableSource);
        panel.add(new JScrollPane(table));
        //trying to get table to redisplay correctly when 
        //changing whichValue.  Still not working correctly.
        panel.invalidate();
        panel.validate();
        doUpdate();
        panel.repaint();        
    }
        
    /**
     * Mutator method for whether table should include a column of "x" values.
     */
    public void setShowXColumn(boolean b) {
        showXColumn = b;
        setupDisplay();
    }
    /**
     * Accessor method for flag indicating if table should include a column of "x" values.
     */
    public boolean isShowXColumn() {return showXColumn;}
    
    public void repaint() {table.repaint();}

    public Component graphic(Object obj) {return panel;}

    private class MyTableData extends AbstractTableModel {
        
        String[] columnNames;
        Class[] columnClasses;
        int nColumns;
        int y0;  //index of first y column (1 if showing x, 0 otherwise)
        
        MyTableData() {
            nColumns = ySource.length;
            y0 = showXColumn ? 1 : 0;
            nColumns += y0;
            columnNames = new String[nColumns];
            columnClasses = new Class[nColumns];
            if(showXColumn) {
                if(xSource != null) {
                    if(xSource instanceof DataSource.X) columnNames[0] = ((DataSource.X)xSource).getXLabel();
                    else columnNames[0] = xSource.getLabel();
                }
                else columnNames[0] = "";
                columnClasses[0] = Double.class;
            }
            for(int i=y0; i<nColumns; i++) {
                columnNames[i] = (ySource[i-y0] != null) ? ySource[i-y0].getLabel() : "";
                columnClasses[i] = Double.class;
            }
        }
        
            //x-y values are updated in update() method
            //because we don't want to call currentValue
            //or average for each function entry
        public Object getValueAt(int row, int column) {
            if(showXColumn && column == 0) return new Double(x[row]);
            else return new Double(y[column-y0][row]);
        }
        
        public int getRowCount() {return (y != null && y[0] != null) ? y[0].length : 0;}
        public int getColumnCount() {return y0 + ((ySource != null) ? ySource.length : 0);}
                
        public String getColumnName(int column) {return columnNames[column];}
        public Class getColumnClass(int column) {return columnClasses[column];}
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        Default.ATOM_SIZE = 1.0;                   
	    IntegratorHard integratorHard = new IntegratorHard();
	    SpeciesDisks speciesDisks = new SpeciesDisks();
	    speciesDisks.setNMolecules(300);
	    Phase phase = new Phase();
	    Potential2 potential = new P2HardSphere();
	    Controller controller = new Controller();
	    DisplayPhase displayPhase = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard.new Timer(integratorHard.chronoMeter());
	    timer.setUpdateInterval(10);
        //part that is unique to this demonstration
        Default.BLOCK_SIZE = 20;
        MeterFunction rdf = new MeterRDF();
        rdf.setActive(true);
        DisplayTableFunction rdfTable = new DisplayTableFunction();
        //end of unique part
        
		Simulation.instance.elementCoordinator.go();
		                                    
        potential.setIterator(new AtomPairIterator(phase));
//        potential.set(species.getAgent(phase));
		
		Simulation.instance.setBackground(java.awt.Color.yellow);
        Simulation.makeAndDisplayFrame(Simulation.instance);
        		                                    
        rdfTable.setDataSources(rdf);
//        rdfTable.setX(rdf.X());
    }//end of main

}