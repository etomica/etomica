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
import jclass.table3.VectorDataSource;
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
public class MeterVolume extends Panel implements JCActionListener, IntegrationIntervalListener
{
    transient Vector updateVector = new Vector(5);
    private final Double value = new Double(0.0);
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
    public  ErrorBlock block;
    public  ErrorBlock blockT;
    transient private Double zero;
    transient int cellCount=0;
    transient int integrationInterval;
    
    
    static int PRESSURE = 0;
    static int CONFIGS = 1;
    static int CURRENT = 2;
    static int AVERAGE = 3;
    static int ERROR = 4;

    JCTextField tempField;

    int symbolSize = 10;
    Color configs_bg = new Color(204, 204, 255);	// Table configs BG
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

    
    
    public MeterVolume()
    {
        block = new ErrorBlock();
        blockT = new ErrorBlock();
        zero = new Double(0.0);
        updateVector.addElement(new Double(-1.0));
        updateVector.addElement(new Double(-1.1));
        updateVector.addElement(new Double(-1.2));
        updateVector.addElement(new Double(-1.3));
        updateVector.addElement(new Double(-1.4));
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
        
        add("West",tblPanel);
 //       ((VectorDataSource)tbl.getDataSource()).deleteColumns(1,4);

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
	private AtomWall leftWall, rightWall, bottomWall, piston;
	private int nAtoms;
	private double scale3, scale2, scale, atomDiameter, simAtomDiameter;
	public void updateAverage(IntegrationIntervalEvent iie) {
	    iieCount++;
	    if(firstCall && (iieCount - updateInterval)== 0) {
    	    firstCall = false;
    	    blockCount++;
    	    Vector species = iie.phase.speciesVector;
    	    leftWall = (AtomWall)((Species)species.elementAt(1)).firstAtom;
    	    bottomWall = (AtomWall)((Species)species.elementAt(2)).firstAtom;
    	    rightWall = (AtomWall)((Species)species.elementAt(3)).firstAtom;
    	    piston = (AtomWall)((Species)species.elementAt(4)).firstAtom;
    	    energy = (bottomWall.r[1]-piston.r[1])
    	                *(rightWall.r[0]-leftWall.r[0]);
            block.sum(energy);
            blockT.sum(iie.phase.getKineticTemperature());
            iieCount = 0;
            nAtoms = iie.phase.nAtomTotal;
            atomDiameter = 3.5;
            simAtomDiameter = ((SpeciesDisk)species.elementAt(0)).diameter;
            scale = atomDiameter/simAtomDiameter;
            scale2 = scale*scale;
            scale3 = scale2*scale;
	    }
	    else if(iieCount == updateInterval) {
	        iieCount = 0;
    	    blockCount++;
    	    energy = nAtoms/(((bottomWall.r[1]-piston.r[1]-simAtomDiameter/2.)*(rightWall.r[0]-leftWall.r[0]-simAtomDiameter/2.)*scale2)*atomDiameter);
            energy *= 1e30/1e3/Constants.AVOGADRO;
            double temperature = iie.phase.getKineticTemperature();
            double pressure = piston.f[1]/(piston.diameter*scale2*atomDiameter)/Constants.BAR2SIM;
            block.sum(energy);
            blockT.sum(temperature);
            if(blockCount % nBlock == 0) {
                block.update();
                blockT.update();
    //            value.doubleValue(pressure);
                updateVector.setElementAt(new Double(pressure),0);
                modelEnergy.setTableDataItem(new Double(pressure), 0, 0);
    //            value.doubleValue(temperature);
                updateVector.setElementAt(new Double(temperature),1);
                modelEnergy.setTableDataItem(new Double(temperature), 1, 0);
    //            value.doubleValue(energy);
                updateVector.setElementAt(new Double(energy),2);
                modelEnergy.setTableDataItem(new Double(energy), 2, 0);
    //            value.doubleValue(blockT.average());
                updateVector.setElementAt(new Double(blockT.average()),3);
                modelEnergy.setTableDataItem(new Double(blockT.average()), 3, 0);
     //           value.doubleValue(block.average());
                updateVector.setElementAt(new Double(block.average()),4);                
                modelEnergy.setTableDataItem(new Double(block.average()), 4, 0);
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
        tbl.setAlignment(JCTblEnum.ALL, PRESSURE, JCTblEnum.BOTTOMCENTER);
        tbl.setBackground(bg_color);
        tbl.setBackground(JCTblEnum.ALL,JCTblEnum.ALL, bg_color);
        tbl.setForeground(JCTblEnum.ALL,JCTblEnum.ALL, data_color);
        tbl.setBackground(JCTblEnum.ALLCELLS, PRESSURE, current_bg);
        tbl.setBackground(JCTblEnum.ALLCELLS, CONFIGS, current_bg);
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
        tbl.setEditable(JCTblEnum.ALLCELLS, CURRENT, false);
        tbl.setEditable(JCTblEnum.ALLCELLS, AVERAGE, false);
        tbl.setEditable(JCTblEnum.ALLCELLS, ERROR, false);
        tbl.setJumpScroll(JCTblEnum.JUMP_VERTICAL);
        tbl.setRowLabelDisplay(false);
        tbl.setCellBorderWidth(1);
        tbl.setTraversable(JCTblEnum.ALL,JCTblEnum.ALL,true);
        tbl.setTraversable(JCTblEnum.ALLCELLS, CURRENT, true);
        tbl.setTraversable(JCTblEnum.ALLCELLS, AVERAGE, true);
        tbl.setTraversable(JCTblEnum.ALLCELLS, ERROR, true);
        tbl.setPixelWidth(JCTblEnum.ALLCELLS, JCTblEnum.VARIABLE);
 //       tbl.setVertSBDisplay(JCTblEnum.SBDISPLAY_ALWAYS);
    }
    
    /**
    * Initialize basic chart look and feel
    */
    public void setupChart(JCChart chart) {
        chart.setAllowUserChanges(false);
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