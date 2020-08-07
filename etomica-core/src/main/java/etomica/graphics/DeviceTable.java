/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.table.TableModel;
import java.awt.*;

/**
 * Presents a table of numeric properties that can be edited by typing in values.
 *
 */
public class DeviceTable /*implements EtomicaElement*/ {
    
    private JTable table;
    private JPanel panel;
    private JScrollPane scrollPane;

    public DeviceTable(String[] columnNames) {
        this(new DeviceTableModelGeneric(null, columnNames));
    }

    public DeviceTable(TableModel tableModel) {
        super();

        panel = new JPanel();
        panel.setBorder(new TitledBorder(null, "Table", TitledBorder.CENTER, TitledBorder.TOP));
        panel.setLayout(new GridBagLayout());
//        panel.setSize(new java.awt.Dimension(350, 100));
//        panel.setPreferredSize(new java.awt.Dimension(350, 100));

        table = new JTable(tableModel);
//        table.setFillsViewportHeight(true);
        table.setSelectionBackground(java.awt.Color.YELLOW);

        initCellEditor(tableModel);

        scrollPane = new JScrollPane(table);

        GridBagLayout gbLayout = new GridBagLayout();
        GridBagConstraints gbConst = new GridBagConstraints();
        gbConst.gridx = 0; gbConst.gridy = 0;
        gbLayout.setConstraints(table, gbConst);
        
//        scrollPane.setSize(new java.awt.Dimension(400, 100));
//        scrollPane.setPreferredSize(new java.awt.Dimension(400, 100));

        panel.add(scrollPane);

    }

    public void initCellEditor(TableModel tableModel) {
	    for(int col = 0; col < tableModel.getColumnCount(); col++) {
	        table.getColumn(tableModel.getColumnName(col)).setCellEditor(new ValueEditor());
	    }
    }

    /**
     * Set the preferred size attribute of the scrolled window the table
     * sits on.
     * @param width
     * @param height
     */
    public void setPreferredSize(int width, int height) {
    	scrollPane.setPreferredSize(new java.awt.Dimension(width, height));
    }

    /**
     * Set the size attribute of the scrolled window the table sits on.
     * @param width
     * @param height
     */
    public void setSize(int width, int height) {
    	scrollPane.setSize(width, height);
    }

    /**
     * Set the title of the border surrounding the table.
     * @param title
     */
    public void setTitle(String title) {
        ((TitledBorder)panel.getBorder()).setTitle(title);
    }

    /**
     * Returns the indices of all the selected rows on the table.
     * @return an array of integers containing the indices of all selected rows,
     * or an empty array if no row is selected
     */
    public int[] getSelectedRows() {
    	return table.getSelectedRows();
    }

    /**
     * Returns the top level panel that the table components sit on.
     */
    public Component graphic() {
        return panel;
    }

    private class ValueEditor extends DefaultCellEditor {
        
        public ValueEditor() {
            super(new javax.swing.JTextField());
            setClickCountToStart(2);
        }

        public boolean isCellEditable(java.util.EventObject e) {
        	return super.isCellEditable(e);
        }

    }//end of ValueEditor


    /**
     * main method to demonstrate and test this class.
     */
    public static void main(String[] args) {
        
        javax.swing.JFrame f = new javax.swing.JFrame();
        f.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
//        f.setPreferredSize(new java.awt.Dimension(600, 350));
        f.setSize(new java.awt.Dimension(600, 350));

        etomica.graphics.DeviceSlider slide1 = new etomica.graphics.DeviceSlider(null);
        etomica.graphics.DeviceSlider slide2 = new etomica.graphics.DeviceSlider(null);

        String p1 = "0.0";
        String p2 = "0.1";
        String p3 = "0.2";
        String p4 = "1.0";
        String p5 = "1.1";
        String p6 = "1.2";
        
        DeviceTableModelGeneric td = new DeviceTableModelGeneric(new Object[][] {{p1, p2, p3}, {p4, p5, p6}}, new String[] {"X", "Y", "Z"});
//      DeviceTableModel td = new DeviceTableModel(null, new String[] {"X", "Y", "Z"});
        DeviceTable table = new DeviceTable(td);
        table.setPreferredSize(400, 75);

        javax.swing.JPanel mainPanel = new javax.swing.JPanel();
        mainPanel.setLayout(new java.awt.GridLayout(3, 1));
        mainPanel.add(slide1.graphic());
        mainPanel.add(slide2.graphic());
        mainPanel.add(table.graphic());
        f.getContentPane().add(mainPanel);

        f.pack();
        f.setVisible(true);
    }//end of main

    
}



