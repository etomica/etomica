package simulate;

import java.awt.*;
import java.util.Vector;
import java.util.*;
import jclass.chart.ChartDataView;
import jclass.chart.JCChartStyle;
import jclass.chart.JCSymbolStyle;
import jclass.chart.EventTrigger;
import jclass.chart.JCChart;
import jclass.chart.FileDataSource;
import jclass.chart.SimpleChart;
import jclass.chart.Chartable;
import jclass.chart.*;

import com.sun.java.swing.*;
public class rdfChart extends java.awt.Panel implements IntegrationIntervalListener
{

	public rdfChart()
	{
		//{{INIT_CONTROLS
		setLayout(null);
		setSize(200,200);
		simpleChart1 = new jclass.chart.SimpleChart();
		simpleChart1.setMargins(new Insets(1,1,1,1));
		simpleChart1.setBounds(130,12,184,183);
		add(simpleChart1);
		button1 = new java.awt.Button();
		button1.setLabel("Reset");
		button1.setBounds(338,44,48,23);
		button1.setBackground(new Color(12632256));
		add(button1);
		//}}
		
		rdf1 = new rdf();
		
	    simpleChart1.setAllowUserChanges(true);
        simpleChart1.setTrigger(0, new EventTrigger(Event.META_MASK,EventTrigger.CUSTOMIZE));
	    rdfView = simpleChart1.getDataView(0);
	    rdfView.setDataSource(rdf1);
	    style=rdfView.getChartStyle(0);
	    style.setSymbolShape(JCSymbolStyle.NONE);
	
//	    initializeRDF();
	
	}
	//{{DECLARE_CONTROLS
	jclass.chart.SimpleChart simpleChart1;
	java.awt.Button button1;
	//}}
	JCChartStyle style;
    ChartDataView rdfView;
    rdf rdf1;
    int updateInterval = 1;
    int iieCount = 0;
    int nSum = 0;
    int nPoints = 100;
    double[] rdfSum;
    double rMin = 0.0, rMax = 0.5, deltaR;
    double[] r, shellVolume;
    double density;
    boolean rMaxHalfEdgeLength = true;  //flag indicating that rMax is computed as half the edgelength of the simulation box
    boolean firstCall = true;
    
    public int getUpdateInterval() {return updateInterval;}
    public void setUpdateInterval(int i) {updateInterval = i;}
    
    public int getNPoints() {return nPoints;}
    public void setNPoints(int i) {
        nPoints = i;
        initializeRDF();
    }
    
    public double getRMin() {return rMin;}
    public void setRMin(double r) {
        rMin = r;
        initializeRDF();
    }
	
    public double getRMax() {return rMax;}
    public void setRMax(double r) {
        rMax = r;
        initializeRDF();
    }
    
    public boolean isRMaxHalfEdgeLength() {return rMaxHalfEdgeLength;}
    public void setRMaxHalfEdgeLength(boolean b) {rMaxHalfEdgeLength = b;}
	
	private void initializeRDF() {
	    deltaR = (rMax - rMin)/(double)nPoints;
	    r = new double[nPoints];
	    rdfSum = new double[nPoints];
	    shellVolume = new double[nPoints];
	    Vector rowR = new Vector();
	    Vector rowRDF = new Vector();
	    nSum = 0;
	    for(int i=0; i<nPoints; i++) {
	        r[i] = ((double)i + 0.5)*deltaR;
	        rdfSum[i] = 0;
            rowR.addElement(new Double(r[i]));
            rowRDF.addElement(new Double(1.0));
	        shellVolume[i] = Math.PI*( (i+1)*(i+1) - i*i )*deltaR*deltaR;
	    }
         rdf1.setDataRow(rowR,0);
         rdf1.notifyObservers(new ChartDataModelUpdate(ChartDataModelUpdate.CHANGE_ROW, 0, -1));
         rdf1.setDataRow(rowRDF,1);
         rdf1.notifyObservers(new ChartDataModelUpdate(ChartDataModelUpdate.CHANGE_ROW, 1, -1));
	}
		
	public void updateAverage(IntegrationIntervalEvent iie) {
	    if(firstCall) {
	        firstCall = false;
	        if(rMaxHalfEdgeLength) {
	            setRMax(0.5);
	        }
	        initializeRDF();
	        density = (double)((Species)iie.phase.speciesVector.elementAt(0)).nElements;
	    }
	    iieCount++;
	    if(iieCount == updateInterval) {
	        iieCount = 0;
	        Phase phase = iie.phase;
	        double[] rdfConfig = calculateDistribution(phase);
	        nSum++;
	        System.out.println("update " + nSum);
	        Vector rowRDF = new Vector();
	        double rdfAverage;
	        for(int i=0; i<nPoints; i++) {
	            rdfSum[i] += rdfConfig[i];
	            rdfAverage = rdfSum[i]/((double)nSum*density*shellVolume[i]);
	            rowRDF.addElement(new Double(rdfAverage));
	        }	            
	        rdf1.setDataRow(rowRDF,1);
            rdf1.notifyObservers(new ChartDataModelUpdate(ChartDataModelUpdate.CHANGE_ROW, 1, -1));
	        
	    }
	}
	 
	public double[] calculateDistribution(Phase phase) {
	    double r2, rMax2 = rMax*rMax;
	    double[] sum = new double[nPoints];
	    int nCenters;
	    nCenters = 0;
	    for(int j=0; j<nPoints; j++) {sum[j] = 0.0;}
	    for(SpeciesElement e1=phase.firstElement; e1!=null; e1=e1.getNext()) {
	        if(e1 instanceof Molecule) {
	            nCenters++;
	            for(SpeciesElement e2=e1.getNext(); e2!=null; e2=e2.getNext()) {
	                if(e2 instanceof Molecule) {
// update for space	                    r2 = phase.space.rr((Molecule)e1,(Molecule)e2);
	                    r2 = 1.0;  //update
	                    if(r2 < rMax2) {
	                        int j = (int)Math.floor(Math.sqrt(r2)/deltaR);
	                        sum[j] += 1.0;
	                    }
	                }
	            }
	        }
	    }
	    for(int j=0; j<nPoints; j++) {sum[j] *= (2.0/(double)nCenters);}
	    return sum;
	}
	
    class rdf extends jclass.chart.ChartDataModel {

        Vector data;

        rdf() {
            data = new Vector(2);
            data.addElement(new Double(0.0));
            data.addElement(new Double(0.0));
        }

       public void setDataRow(Vector row, int indx) {
         data.setElementAt(row,indx);
         setChanged();
       }
       
       public Object getDataItem(int row, int column) {
	        Object rval = null;
	        try {
		        rval = ((Vector)data.elementAt(row)).elementAt(column);
	        }
	        catch (Exception e) {
	        }
	        return rval;
        }

        public boolean setDataItem(int row, int column, java.lang.Object item) {return false;}


        /**
        * Overridden from Chartable
        */
        public Vector getRow(int row) {
        	Vector rval = null;
        	try {
        		rval = (Vector)data.elementAt(row);
        	}
        	catch (Exception e) {
        	}
        	return rval;
        }       

        /**
         * Overridden from Chartable
         */
        public int getDataInterpretation() {
	        return GENERAL;
        }

        /**
         * Overridden from Chartable
        */
        public int getNumRows() {
	        if (data == null) return 0;
	        return data.size();
        }

        /**
        * Overridden from Chartable
        */
        public String[] getPointLabels() {
	        return null;
        }

        /**
         * Overridden from Chartable
         */
        public String getSeriesName(int row) {
        	return null;
        }

        /**
         * Overridden from Chartable
         */
        public String getSeriesLabel(int row) {
        	return getSeriesName(row);
        }

        /**
        * Overridden from Chartable
        */
        public String getName() {
	        return new String("Dave Data");
        }
        
    }  //end of class rdf

	static public void main(String args[])
	{
		class DriverFrame extends java.awt.Frame {
			public DriverFrame() {
				addWindowListener(new java.awt.event.WindowAdapter() {
					public void windowClosing(java.awt.event.WindowEvent event)
					{
						dispose();	  // free the system resources
						System.exit(0); // close the application
					}
				});
				this.setLayout(new java.awt.BorderLayout());
				this.setSize(300,300);
				this.add(new rdfChart());
			}
		}

		new DriverFrame().show();
	}


}