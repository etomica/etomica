package simulate;

import simulate.*;
import java.beans.*;
import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.JTable;
import com.sun.java.swing.table.*;
import com.sun.java.swing.Box;

public class DisplayTable extends simulate.Display
{
    public JTable table;
    MyTableData dataSource;
    Meter[] meter = new Meter[0];
    int nMeters = 0;
    private boolean showingPE = false;  //flag to indicate if PotentialEnergy meter should be displayed
    private boolean showingKE = false;  //same for KineticEnergy meter
    Box panel = Box.createVerticalBox();
    Button resetButton = new Button("Reset averages");
    
    public DisplayTable()
    {
        super();
        dataSource = new MyTableData();   //inner class, defined below
        table = new JTable(dataSource);
		resetButton.addActionListener(new ActionListener() {
		    public void actionPerformed(java.awt.event.ActionEvent event) {DisplayTable.this.resetAverages();}
		});
//        setVisible(false);  //get out of way of painting of table
        panel.add(table);
        panel.add(resetButton);
        add(panel);
    }
    
    public void resetAverages() {
        for(int i=0; i<nMeters; i++) {meter[i].reset();}
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
        for(int i=0; i<meter.length; i++) {temp[i] = meter[i];}
        temp[nMeters-1] = m;
        meter = temp;
    }
        
    public void setPhase(Phase p) {
        for(Meter m=p.firstMeter; m!=null; m=m.nextMeter()) {this.addMeter(m);}
    }

    public void doUpdate() {}
    public void repaint() {table.repaint();}
    
    class MyTableData extends AbstractTableModel {
        
        String[] columnNames = new String[] {"Property", "Current", "Average"};
        Class[] columnClasses = {String.class, Double.class, Double.class};
        
        MyTableData() {}
        
        public Object getValueAt(int row, int column) {
            Meter m = meter[row];
            switch(column) {
                case 0: return m.getLabel();
                case 1: return new Double(m.currentValue());
                case 2: return new Double(m.average());
                default: return null;
            }
        }
        
        public int getRowCount() {return meter.length;}
        public int getColumnCount() {return 3;}  //label, current, average
                
        public String getColumnName(int column) {return columnNames[column];}
        public Class getColumnClass(int column) {return columnClasses[column];}
    }
    
}