/*
 * (c) Copyright 1996, 1997, 1998 KL GROUP INC.  
 * ALL RIGHTS RESERVED 
 *
 * This file is provided for demonstration and educational uses only.
 * Permission to use, copy, modify and distribute this file for
 * any purpose and without fee is hereby granted, provided that the
 * above copyright notice and this permission notice appear in all
 * copies, and that the name of KL Group not be used in advertising
 * or publicity pertaining to this material without the specific,
 * prior written permission of an authorized representative of
 * KL Group.
 *
 * KL GROUP MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES, EX-
 * PRESS OR IMPLIED, WITH RESPECT TO THE SOFTWARE, INCLUDING, BUT
 * NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR ANY PARTICULAR PURPOSE, AND THE WARRANTY AGAINST IN-
 * FRINGEMENT OF PATENTS OR OTHER INTELLECTUAL PROPERTY RIGHTS.  THE
 * SOFTWARE IS PROVIDED "AS IS", AND IN NO EVENT SHALL KL GROUP OR
 * ANY OF ITS AFFILIATES BE LIABLE FOR ANY DAMAGES, INCLUDING ANY
 * LOST PROFITS OR OTHER INCIDENTAL OR CONSEQUENTIAL DAMAGES 
 * RELATING TO THE USE OF THIS SOFTWARE.
 */

//   RCSID -- $RCSfile$ $Revision$
//            $Date$  $Locker$  KL Group Inc.

package simulate;

import java.awt.*;
import java.util.Vector;
import jclass.chart.JCAxis;
import jclass.base.Border;
import jclass.bwt.JCActionEvent;
import jclass.bwt.JCActionListener;
import jclass.bwt.JCGridLayout;
import jclass.bwt.JCLabel;
import jclass.bwt.JCTextField;
import jclass.chart.EventTrigger;
import jclass.chart.JCAxis;
import jclass.chart.JCChart;
import jclass.chart.JCChartStyle;
import jclass.chart.JCLegend;
import jclass.table3.JCTblEnum;
import jclass.table3.StandardCellBorder;
import jclass.table3.Table;
import jclass.util.JCString;

/**
 * A Class which calculates and displays the average and error of the internal energy.
 * This class should be registered as an IntegrationIntervalListener.  When this event
 * is fired from the integrator, this class will increment a counter.  When the counter
 * is equal to a preset value, this class retrives the instaneous energy from the integrator
 * an includes the value in a running average and error.
 *
 * @author David A. Kofke
 * @version %I%, %G%
 */
public class MeterEnergy extends Panel implements JCActionListener, IntegrationIntervalListener
{
    transient Vector updateVector = new Vector(4);
    Table tbl = new Table();
    EnergyDataModel modelEnergy = new EnergyDataModel();
    EnergyTableCells energyTableCells = new EnergyTableCells(modelEnergy);
    JCChart chart = new JCChart();

    int updateInterval = 1;
    int iieCount = 0;
    int blockCount = 0;
    int nBlock = 50;
    boolean firstCall = true;
    transient double initialEnergy;
    transient double energy;
    transient ErrorBlock block;
    transient private Double zero;
    transient int cellCount=0;
    transient int integrationInterval;
    
    
    
    static int CONFIGS = 0;
    static int CURRENT = 1;
    static int AVERAGE = 2;
    static int ERROR = 3;

    JCTextField tempField;

    int symbolSize = 10;
    Color configs_bg = new Color(86, 196, 124);	// Table configs BG
    Color configs_fg = new Color(0, 60, 30);         // Table configs FG
    Color current_bg = new Color(204, 204, 255);	// Table current BG, Light, pale purple
    Color current_fg = new Color(72, 72, 108);        // Table current FG
    Color average_bg = new Color(184, 184, 235);	// Table average BG
    Color average_fg = new Color(0, 60, 30);         // Table average FG
    Color error_bg = new Color(204, 204, 255);	// Table error BG, Light, pale purple
    Color error_fg = new Color(0, 0, 255);        // Table errors FG
	//Color data_color = new Color(0xCC,0xCC,0xFF);	// Nice Light Blue
	Color data_color = new Color(0,0,0);
	
    Color bg_color = new Color(0xD3,0xD3,0xD3);	// Light Gray
    Color grid_color = Color.gray;
    Font  headline_font = new Font("TimesRoman", Font.BOLD, 18);

    
    
    public MeterEnergy()
    {
        block = new ErrorBlock();
        zero = new Double(0.0);
        updateVector.addElement(new Double(-1.0));
        updateVector.addElement(new Double(-1.1));
        updateVector.addElement(new Double(-1.2));
        updateVector.addElement(new Double(-1.3));
	    setBackground(bg_color);
	    BorderLayout layout = new BorderLayout(10,10);
        //GridBagLayout layout = new GridBagLayout();
	    super.setLayout(layout);
        //GridBagConstraints constraints = new GridBagConstraints();
        //constraints.weightx = 100;
        //constraints.weighty = 100;
    
        JCString headline = JCString.parse(this, "[FONT=TimesRoman-bold-16]Simulation Data");
        JCLabel title = new JCLabel(headline);
        add("North", title);

	     // Set up table 
        tbl.setDataSource(energyTableCells);
        this.setupTable(tbl);
        Panel tblPanel = new Panel();
        tblPanel.setLayout(new BorderLayout());
        tblPanel.add("Center",tbl);
        
        // Set up chart
        JCChart chart = new JCChart(JCChart.PLOT);
        //above constructor executes line below automatically.
        //chart.getDataView(0).setChartType(JCChart.PLOT);
        chart.getDataView(0).setDataSource(new EnergyChartView(modelEnergy));
        //chart.getDataView(1).setDataSource(new ErrorChartView(modelEnergy));
        this.setupEnergyChart(chart);

        Panel chartPanel = new Panel();
        //chartPanel.setLayout(new GridLayout(2, 1));
        chartPanel.setLayout(new GridLayout(1, 1));
        chartPanel.add(chart);
        /*
        LabelLast label = new LabelLast(modelEnergy);
        label.setAttachMethod(JCChartLabel.ATTACH_DATAINDEX);
        label.setDataIndex(new JCDataIndex(chart.getDataView(0), chart.getDataView(0).getSeries(0), 0, chart.getDataView(0).getSeries(0).getLastPoint()));
        label.setAnchor(JCChartLabel.EAST);
        label.setOffset(new java.awt.Point(2, 0));
        label.setBorderType(Border.ETCHED_OUT);
        label.setBorderWidth(2);
        chart.addChartLabel(label);

        // Use a custom label generator to format the y axis labels
        chart.getDataView(0).getYAxis().setLabelGenerator(label);
        */        

        Panel tempPanel = new Panel();
        tempPanel.setLayout(new JCGridLayout(1, 2));
        tempField = new JCTextField("300.0", 5);
        tempField.addActionListener(this);
        tempField.setFont(headline_font);
        tempField.setBackground(Color.white);
        tempPanel.add(tempField);
        JCLabel tempLabel = new JCLabel("Temperature (K)");
        tempLabel.setFont(headline_font);
        tempPanel.add(tempLabel);

        Panel topPanel = new Panel();
        topPanel.setLayout(new JCGridLayout(1,2));
        topPanel.add(new JCLabel("      "));
        topPanel.add(tempPanel);
        tblPanel.add("North", topPanel);

        JCString note = JCString.parse(this, "[FONT=TimesRoman-italic-10]Model by [COLOR=blue][HREF=http://www.acsu.edu/~cdbarnes]Dan Barnes[/HREF][/COLOR]");
        JCLabel foot = new JCLabel(note);
        Panel footnote = new Panel();
        footnote.setLayout(new JCGridLayout(1,2));
        footnote.add(foot);
        tblPanel.add("South", footnote);

        add("West",tblPanel);
        add(chartPanel);

        // Trigger initialization
        //this.actionPerformed(null);        
    }
/*    
    public void add(Component c, GridBagConstraints gbc, int x, int y, int w, int h) {
        gbc.gridx = x;
        gbc.gridy = y;
        gbc.gridwidth = w;
        gbc.gridheight = h;
        add(c, gbc);
        System.out.println("x= " + gbc.gridx + " y = " + gbc.gridy + " w = " + gbc.gridwidth + " h = " + gbc.gridheight);
    }
*/    
        
   /**
     * Method returns the number of times the event must fire
     * before the energy is added to the running sum.
     *
     * @return Interval of events to occur before sum is updated with new energy.
     */
    public int getUpdateInterval() {return updateInterval;}
   /**
     * Method sets the number of times the event must fire
     * before the energy is added to the running sum.
     *
     * @parm Interval of events to occur before sum is updated with new energy.
     */
    public void setUpdateInterval(int i) {updateInterval = i;}

   /**
     * This method returns the number of data values a block contains.
     *
     * @parm number of values in one block
     * @see ErrorBlock
     */
    public int getNBlock() {return nBlock;}

   /**
     * This method set the number of data values a block contains.
     *
     * @parm number of values in one block
     * @see ErrorBlock
     */
    public void setNBlock(int nBlock) {this.nBlock = nBlock;}

    public int getSymbolSize() {return symbolSize;}
    
    public void setSymbolSize(int s) {
        symbolSize = s;
        this.setupEnergyChart(chart);
    }
    
   /**
     * This method updates all stored values.  The averages and errors
     * are stored in a Vector is limited to a maximum size of 500.
     *
     * @parm IntegrationIntervalEvent from the Intgrator class
     * @see Integrator
     */
	public void updateAverage(IntegrationIntervalEvent iie) {
	    iieCount++;
	    if(firstCall && (iieCount - updateInterval)== 0) {
    	    firstCall = false;
    	    blockCount++;
    	    energy = iie.phase.getKineticEnergy();
            block.sum(energy);
            iieCount = 0;
	    }
	    else if(iieCount == updateInterval) {
	        iieCount = 0;
    	    blockCount++;
    	    energy = iie.phase.getKineticEnergy();
            block.sum(energy);
            if(blockCount % nBlock == 0) {
                block.update();
                updateVector.setElementAt(new Double(cellCount*(updateInterval*nBlock*(((Integrator)iie.getSource()).getIntegrationInterval()))),0);
                updateVector.setElementAt(new Double(energy),1);
                updateVector.setElementAt(new Double(block.average()),2);
                updateVector.setElementAt(new Double(block.error()),3);
                if(cellCount < 5) {
                    modelEnergy.setTableDataItem(new Double(cellCount*(updateInterval*nBlock*(((Integrator)iie.getSource()).getIntegrationInterval()))), 0, cellCount);
                    modelEnergy.setTableDataItem(new Double(energy), 1, cellCount);
                    modelEnergy.setTableDataItem(new Double(block.average()), 2, cellCount);
                    modelEnergy.setTableDataItem(new Double(block.error()), 3, cellCount);
                }
                else {
                    modelEnergy.addColumn(JCTblEnum.MAXINT, null, updateVector);
                }
                cellCount++;
            }
        }
	}
    //method to resolve VCafe problem of instantiating a new layout manager in the Applet at design time.
    //This forces the original instance to be used (dummy is that generated by the Applet).
    public void setLayout(LayoutManager dummy) {}

    /**
    * Initialize table look and feel
    */
    public void setupTable(Table tbl) {
        tbl.setAutoScroll(JCTblEnum.AUTO_SCROLL_BOTH);
        tbl.setAlignment(JCTblEnum.ALL,JCTblEnum.ALL, JCTblEnum.MIDDLERIGHT);
        tbl.setAlignment(JCTblEnum.LABEL, JCTblEnum.ALL, JCTblEnum.BOTTOMCENTER);
        tbl.setAlignment(JCTblEnum.ALL, CONFIGS, JCTblEnum.BOTTOMCENTER);
        tbl.setBackground(bg_color);
        tbl.setBackground(JCTblEnum.ALL,JCTblEnum.ALL, bg_color);
        tbl.setForeground(JCTblEnum.ALL,JCTblEnum.ALL, data_color);
        tbl.setBackground(JCTblEnum.ALLCELLS, CURRENT, current_bg);
        tbl.setBackground(JCTblEnum.ALLCELLS, AVERAGE, average_bg);
        tbl.setBackground(JCTblEnum.ALLCELLS, ERROR, error_bg);
        tbl.setCellBorderType(JCTblEnum.ALLCELLS, JCTblEnum.ALLCELLS, new StandardCellBorder(JCTblEnum.BORDER_NONE));
        tbl.setCellBorderType(JCTblEnum.LABEL, JCTblEnum.ALL, new StandardCellBorder(JCTblEnum.BORDER_OUT));
        tbl.setCellBorderType(JCTblEnum.ALLCELLS, CURRENT, new StandardCellBorder(JCTblEnum.BORDER_IN));
        tbl.setCellBorderType(JCTblEnum.ALLCELLS, AVERAGE, new StandardCellBorder(JCTblEnum.BORDER_IN));
        tbl.setCellBorderType(JCTblEnum.ALLCELLS, ERROR, new StandardCellBorder(JCTblEnum.BORDER_IN));
        tbl.setCharHeight(JCTblEnum.LABEL, 2);
        tbl.setCharWidth(JCTblEnum.ALL, 5);
        tbl.setFont(JCTblEnum.LABEL, JCTblEnum.ALL, new Font("TimesRoman", Font.BOLD, 12));
        tbl.setFrameBorderType(JCTblEnum.BORDER_OUT);
        tbl.setFrameBorderWidth(1);
        tbl.setEditable(JCTblEnum.ALL,JCTblEnum.ALL,false);
        tbl.setEditable(JCTblEnum.ALLCELLS, CURRENT, true);
        tbl.setEditable(JCTblEnum.ALLCELLS, AVERAGE, true);
        tbl.setEditable(JCTblEnum.ALLCELLS, ERROR, true);
        tbl.setJumpScroll(JCTblEnum.JUMP_VERTICAL);
        tbl.setRowLabelDisplay(false);
        tbl.setCellBorderWidth(1);
        tbl.setTraversable(JCTblEnum.ALL,JCTblEnum.ALL,false);
        tbl.setTraversable(JCTblEnum.ALLCELLS, CURRENT, true);
        tbl.setTraversable(JCTblEnum.ALLCELLS, AVERAGE, true);
        tbl.setTraversable(JCTblEnum.ALLCELLS, ERROR, true);
        tbl.setVertSBDisplay(JCTblEnum.SBDISPLAY_ALWAYS);
    }
    
    /**
    * Initialize basic chart look and feel
    */
    public void setupChart(JCChart chart) {
        chart.setAllowUserChanges(true);
        chart.setBackground(bg_color);
        chart.getChartArea().getPlotArea().setBackground(Color.lightGray);
        chart.getChartArea().setAxisBoundingBox(true);
        chart.getChartArea().setFastAction(true);
        chart.getHeader().setIsShowing(true);
        chart.getHeader().setTop(0);

        chart.getLegend().setIsShowing(true);
        chart.getLegend().setAnchor(JCLegend.NORTHWEST);
        chart.getLegend().setBorderType(Border.ETCHED_OUT);

        JCAxis xaxis = chart.getDataView(0).getXAxis();
        xaxis.setAnnotationMethod(JCAxis.POINT_LABELS);
        xaxis.setGridIsShowing(true);
        xaxis.getTitle().setText("[FONT=Times-bold-14]Configurations", true);
        xaxis.getGridStyle().getLineStyle().setColor(grid_color);

        JCAxis yaxis = chart.getDataView(0).getYAxis();
        yaxis.setPlacement(JCAxis.MIN);
        yaxis.setGridIsShowing(true);
        yaxis.getGridStyle().getLineStyle().setColor(grid_color);

        /*
        JCSymbolStyle.NONE, JCSymbolStyle.DOT,
        JCSymbolStyle.BOX, JCSymbolStyle.TRIANGLE,
        JCSymbolStyle.DIAMOND, JCSymbolStyle.STAR,
        JCSymbolStyle.VERT_LINE, JCSymbolStyle.HORIZ_LINE
        JCSymbolStyle.CROSS, JCSymbolStyle.CIRCLE
        JCSymbolStyle.SQUARE
        */
        JCChartStyle style = new JCChartStyle();
        style.setSymbolShape(jclass.chart.JCSymbolStyle.TRIANGLE);
        style.setSymbolColor(current_fg);
        style.setSymbolSize(symbolSize);
        style.getLineStyle().setColor(current_fg);
        style.setLineWidth(2);
        style.setFillColor(current_fg);
        chart.getDataView(0).getSeries(0).setStyle(style);

        style = new JCChartStyle();
        style.setSymbolShape(jclass.chart.JCSymbolStyle.CIRCLE);
        style.setSymbolColor(average_fg);
        style.setSymbolSize(symbolSize);
        style.getLineStyle().setColor(average_fg);
        style.setLineWidth(2);
        style.setFillColor(average_fg);
        chart.getDataView(0).getSeries(1).setStyle(style);

        style = new JCChartStyle();
        style.setSymbolShape(jclass.chart.JCSymbolStyle.STAR);
        style.setSymbolColor(error_fg);
        style.setSymbolSize(symbolSize);
        style.getLineStyle().setColor(error_fg);
        style.setLineWidth(2);
        style.setFillColor(error_fg);
        chart.getDataView(0).getSeries(2).setStyle(style);
    }


    /**
    * Initialize Energy chart look and feel
    */
    public void setupEnergyChart(JCChart chart) {
        this.setupChart(chart);

        chart.getHeader().setText("Energy Plot[FONT=Times-plain-12][NEWLINE]Drag points", true);
        chart.getHeader().setFont(headline_font);
        chart.getDataView(0).getXAxis().setPlacement(JCAxis.MIN);
        chart.getDataView(0).getYAxis().setIsEditable(false);

        chart.getChartArea().setDepth(3);
        chart.getChartArea().setElevation(12);
        chart.getChartArea().setRotation(20);

        chart.setTrigger(0, new EventTrigger(0, EventTrigger.EDIT));
        chart.setTrigger(0, new EventTrigger(Event.CTRL_MASK, EventTrigger.ROTATE));
        chart.setTrigger(0, new EventTrigger(Event.SHIFT_MASK, EventTrigger.TRANSLATE));
        chart.setTrigger(0, new EventTrigger(Event.META_MASK, EventTrigger.CUSTOMIZE));
    }
    
    public void actionPerformed(JCActionEvent ev) {
        String test = tempField.getText();
        Double value = new Double(test);
        if (value != null) {
            //pass something to simulation.
        }
        else {
           //pass something else to simulation.
        }
    }
}