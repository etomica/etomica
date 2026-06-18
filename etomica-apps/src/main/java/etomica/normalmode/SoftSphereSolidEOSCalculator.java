/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import javax.swing.*;
import java.awt.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;

/**
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class SoftSphereSolidEOSCalculator extends JPanel implements PropertyChangeListener {

	public SoftSphereSolidEOSCalculator() {
		super(new BorderLayout());
		
		double density = 1.1964;
		double temperature = 1.0;
		double[] quantity = getAllQuantity(FCC, n, temperature, density);
		
		//Create the labels.
		densityLabel = new JLabel(densityString);
		temperatureLabel = new JLabel(temperatureString);
		bALabel = new JLabel(bAString);
		ZLabel = new JLabel(ZString);
		bULabel = new JLabel(bUString);
		
		//Create the text fields and set them up.
		densityField = new JFormattedTextField(densityFormat);
		densityField.setValue(density);
		densityField.setColumns(10);
		densityField.addPropertyChangeListener("value", this);
		
		temperatureField = new JFormattedTextField(temperatureFormat);
		temperatureField.setValue(temperature);
		temperatureField.setColumns(10);
		temperatureField.addPropertyChangeListener("value", this);
		
		bAField = new JFormattedTextField(bAFormat);
		bAField.setValue(quantity[0]);
		bAField.setColumns(10);
		bAField.setEditable(false);
//		paymentField.setForeground(Color.red);
		
		ZField = new JFormattedTextField(ZFormat);
		ZField.setValue(quantity[1]);
		ZField.setColumns(10);
		ZField.setEditable(false);
		
		bUField = new JFormattedTextField(bUFormat);
		bUField.setValue(quantity[2]);
		bUField.setColumns(10);
		bUField.setEditable(false);
		
		//Tell accessibility tools about label/textfield pairs.
		densityLabel.setLabelFor(densityField);
		temperatureLabel.setLabelFor(temperatureField);
		bALabel.setLabelFor(bAField);
		ZLabel.setLabelFor(ZField);
		bULabel.setLabelFor(bUField);
		
		//Lay out the labels in a panel.
		JPanel labelPane = new JPanel(new GridLayout(0,1));
		labelPane.add(densityLabel);
		labelPane.add(temperatureLabel);
		labelPane.add(bALabel);
		labelPane.add(ZLabel);
		labelPane.add(bULabel);
		
		//Layout the text fields in a panel.
		JPanel fieldPane = new JPanel(new GridLayout(0,1));
		fieldPane.add(densityField);
		fieldPane.add(temperatureField);
		fieldPane.add(bAField);
		fieldPane.add(ZField);
		fieldPane.add(bUField);
		
		//Put the panels in this panel, labels on left,
		//text fields on right.
		setBorder(BorderFactory.createEmptyBorder(40, 40, 40, 40));
		add(labelPane, BorderLayout.CENTER);
		add(fieldPane, BorderLayout.LINE_END);
		

//		LayoutStyle layout = new 
	}
	
	/** Called when a field's "value" property changes. */
	public void propertyChange(PropertyChangeEvent e) {
		Object source = e.getSource();
		if (source == densityField) {
			density = ((Number)densityField.getValue()).doubleValue();
		} else if (source == temperatureField) {
			temperature = ((Number)temperatureField.getValue()).doubleValue();
		} 
	
		double[] quantity = getAllQuantity(FCC, n, temperature, density);
		bAField.setValue(quantity[0]);
		ZField.setValue(quantity[1]);
		bUField.setValue(quantity[2]);
	}
	
	/**
	* Create the GUI and show it.  For thread safety,
	* this method should be invoked from the
	* event dispatch thread.
	*/
	private static void createAndShowGUI() {
		//Create and set up the window.
		JFrame frame = new JFrame("Soft Sphere EOS Calculator");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		//Add contents to the window.
		frame.add(new SoftSphereSolidEOSCalculator());
		
		//Display the window.
		frame.pack();
		frame.setVisible(true);
	}
	
	public static void main(String[] args) {
		//Schedule a job for the event dispatch thread:
		//creating and showing this application's GUI.
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
			//Turn off metal's use of bold fonts
				UIManager.put("swing.boldMetal", Boolean.FALSE);
				createAndShowGUI();
			}
		});
	}
	
//	public static class Applet extends javax.swing.JApplet{
//		
//		public void init(){
////			createAndShowGUI();
//			getContentPane().add(this);
//		}
//
//		private static final long serialVersionUID = 1L;
//
//	}
	
	public double getUlat0(int id, int n){
		return ulat0[id][n/3-2];
	}
	
	public double getVib0(int id, int n){
		return vib0[id][n/3-2];
	}
	
	public double[] getCoeffAc(int id , int n){
		return c[id][n/3-2];
	}
	
	public double[] getAcQuantity(int id, int n, double temp, double rho){
		double bAc = 0.0;
		double dAcdrho = 0.0;
		double dbAcdbeta = 0.0;
		
		double[] c = getCoeffAc(id, n);
		double upsilon = getUpsilon(n, temp, rho);
		double ups = 1; 
		
		for(int i=1; i<=7; i++){
			ups *= upsilon;
			bAc += c[i]*ups;
			
			dAcdrho += temp*c[i]*ups*(-n*i/(3.0*rho));
			dbAcdbeta += c[i]*ups*(-i*temp);
		}
		
		double[] values = new double[]{bAc, dAcdrho, dbAcdbeta};
		
		return values;
	}
	
	public double[] getAllQuantity(int id, int n, double temp, double rho){
		double upsilon = getUpsilon(n, temp, rho);
		double ulat0 = getUlat0(id, n);
		double Vib0 = getVib0(id, n);
		double[] AcValues = getAcQuantity(id, n, temp, rho); 
		
		double bUlat = ulat0/upsilon; 			            // Equation (A2)
		double bAVib = Vib0 - 3.0/2.0*Math.log(upsilon) 
		                + Math.log(rho*sigma*sigma*sigma);  // Equation (A3)
		double bAc   = AcValues[0];			 				// Equation ( 9)
		double bA    = bUlat + bAVib + bAc; 				// Equation (A1) 
		
		double p = rho*rho*((n/3.0)*(temp/(upsilon*rho))*ulat0  	
				   + AcValues[1]
				   + (n*temp)/(2*rho) + 1/rho);
		
		double u = temp/upsilon*ulat0
				   + AcValues[1]
				   + (3.0/2.0)*temp;
		double[] quantity = new double[]{bA, p/(rho*temp), u/temp};
		return quantity;
	}
	
	public double getUpsilon(int n, double temp, double rho){
		return (temp/epsilon)*Math.pow((rho*sigma*sigma*sigma), -n/3.0);
	}
	
	private static final long serialVersionUID = 1L;
	//Values for the fields
	private double density = 1.1964;
	private double temperature = 1.0;
	private int n = 12;
	
	//Labels to identify the fields
	private JLabel densityLabel;
	private JLabel temperatureLabel;
	private JLabel bALabel;
	private JLabel ZLabel;
	private JLabel bULabel;
	
	//Strings for the labels
	private static String densityString = "Density: ";
	private static String temperatureString = "Temperature: ";
	private static String bAString = "betaA: ";
	private static String ZString = "Z: ";
	private static String bUString = "betaU: ";
	
	//Fields for data entry
	private JFormattedTextField densityField;
	private JFormattedTextField temperatureField;
	private JFormattedTextField bAField;
	private JFormattedTextField ZField;
	private JFormattedTextField bUField;
	
	private NumberFormat densityFormat;
	private NumberFormat temperatureFormat;
	private NumberFormat bAFormat;
	private NumberFormat ZFormat;
	private NumberFormat bUFormat;

    public final static int FCC = 0;
    public final static int HCP = 1;
    public final static int BCC = 2;
    protected double epsilon = 1.0;
    protected double sigma = 1.0;
    
    protected final double[][] ulat0 = new double[][]{{3.613475382,2.208391122,1.516485025},{3.613718378,2.208528128,1.516536721},{3.630711328}};
    protected final double[][] vib0 = new double[][]{{2.69819834,3.48217524,3.88845245},{2.70269996,3.48590058,3.89158992},{2.56473110}};
    protected final double[][][] c = new double[][][]{
    		{{0.0, -0.024270, -0.385254, -1.989540, 24.366184, -208.856577, 860.698810,-1469.431853},
    		 {0.0,  0.249526, -0.245724, -0.500979,  4.258323,  -16.542027,  30.811960,  -23.208566},
    		 {0.0,  0.464148, -0.423821,  0.466322, -0.526094,    0.148125,   0.401893,   -0.384095}},
    		
    		{{0.0, -0.023687, -0.351963, -2.934478, 37.070278, -280.028273, 953.287096, -1140.820226},
    		 {0.0, 0.250099,  -0.282838,  0.007132,  0.922672,   -5.218832,  10.890988,    -7.369508},
    		 {0.0, 0.462065,  -0.413562,  0.429175, -0.461716,    0.111491,   0.368595,    -0.352244}},
    		
    		{{0.0, 0.283913, -2.421458, 14.558777, -77.105689, 229.711416, -391.161702, 335.089915}},
    			 
    };

}
