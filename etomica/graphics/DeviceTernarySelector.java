package etomica.graphics;
import etomica.*;

import java.awt.*;

/**
 * Selector for fractions in a ternary system.  Presents a triangle that can
 * be clicked upon to select (for example) the three mole fractions of a ternary
 * mixture.
 *
 * @author David Kofke
 */
public class DeviceTernarySelector extends Device implements EtomicaElement {
    
    private Triangle triangle;
    private int sideLength, h;
    private int precision;
    private int power; //10^precision
    private String[] labels = new String[3];
    private String[] symbols = labels;
    private boolean showValues = false;
    private javax.swing.JTextField[] boxes;
    private javax.swing.JPanel boxPanel;
    private javax.swing.JPanel trianglePanel;
    private javax.swing.JPanel panel;
    private final SimulationEventManager listenerManager = new SimulationEventManager();
    
    public DeviceTernarySelector() {
        this(new String[] {"A", "B", "C"});
    }
    /**
     * Creates a new devices using the given labels for the diagram.
     * The given string array must be of dimension 3.
     */
    public DeviceTernarySelector(String[] labels) {
        trianglePanel = new TrianglePanel();
        panel = new javax.swing.JPanel();
        panel.add(trianglePanel);
        setLabels(labels);
        setSideLength(200);
        trianglePanel.addMouseListener(new ClickListener());
        setPrecision(3);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Select fractions in a ternary system");
        return info;
    }
    
    public void addListener(DeviceTernarySelector.Listener listener) {
        listenerManager.addListener(listener);
    }
    public void removeListener(DeviceTernarySelector.Listener listener) {
        listenerManager.removeListener(listener);
    }
    
    private int x0() {return panel.getFontMetrics(panel.getFont()).stringWidth(symbols[0]);}
    private int y0() {return panel.getFontMetrics(panel.getFont()).getHeight();}
    
    /**
     * Mutatator method for the precision of the mole fractions.  Precision
     * is the number of significant figures to which the values are rounded.
     * Minimum is 1, maximum is 10.  Default is 2.
     */
    public void setPrecision(int n) {
        n = (n < 1) ? 1 : n; //minimum value is 1
        n = (n > 10) ? 10 : n; //maximum value is 10
        precision = n;
        power = 1;
        for(int i=0; i<n; i++) power *= 10;
    }
    /**
     * Accessor method for the precision of the mole fractions.
     */
    public int getPrecision() {
        return precision;
    }
    
    /**
     * Mutator method for flag indicating whether a set of boxes is
     * displayed with labels and fraction values as text.
     */
    public void setShowValues(boolean b) {
        if(b == showValues) return;//do nothing if not changing value
        showValues = b;
        if(showValues) {
            boxes = new javax.swing.JTextField[3];
            boxPanel = new javax.swing.JPanel();
            boxPanel.setLayout(new GridLayout(3,2));
            panel.add(boxPanel);
            for(int i=0; i<3; i++) {
                boxes[i] = new javax.swing.JTextField(precision+2);
                boxes[i].setEditable(false);
                boxPanel.add(new javax.swing.JLabel(labels[i]));
                boxPanel.add(boxes[i]);
            }
        } else {//turning off boxes
            panel.remove(boxPanel);
            boxes = null;
        }
    }
    /**
     * Accessor method for flag indicating if values are shown in accompanying text boxes.
     */
    public boolean isShowValues() {return showValues;}
    
    private void makeTriangle() {
        triangle = new Triangle(sideLength, x0(), y0());
    }
    
    /**
     * Mutator method for the length of each side of the triangle, in pixels.
     * Default value is 200.
     */
    public void setSideLength(int L) {
        sideLength = L;
        makeTriangle();
    }
    /**
     * Accessor method for the length of each side of the triangle, in pixels.
     * Default value is 200.
     */
    public int getSideLength() {return sideLength;}
    
    /**
     * Labels used to indicate the three component fractions.  These strings
     * are used to label the numeric displays of the fractions, if the device
     * is set to show them.  They are also used to label each vertex of the
     * triangle, if the symbols are not set explicitly.  Default values of
     * the labels are "A", "B", and "C"
     */
    public void setLabels(String[] newLabels) {
        if(newLabels.length != 3) {
            throw new IllegalArgumentException("Invalid dimension in DeviceTernarySelector.setLabels");
        }
        for(int i=0; i<3; i++) {labels[i] = newLabels[i];}
        if(labels == symbols) makeTriangle();
    }
    /**
     * Accessor method for the labels.
     */
    public String[] getLabels() {return labels;}
    
    /**
     * Symbols used to label each vertex of the triangle.  By default the symbols are
     * the same as the labels, and will change any time the labels are changed, unless
     * having once been set explicitly using this method.
     */
    public void setSymbols(String[] newSymbols) {
        if(newSymbols.length != 3) {
            throw new IllegalArgumentException("Invalid dimension in DeviceTernarySelector.setSymbols");
        }
        symbols = new String[3];
        for(int i=0; i<3; i++) {symbols[i] = newSymbols[i];}
        makeTriangle();
    }
    /**
     * Accessor method for the symbols.
     */
    public String[] getSymbols() {return labels;}
    
    /**
     * Graphical display given by this device.
     */
    public java.awt.Component graphic(Object obj) {
        return panel;
    }
    
    /**
     * Simple panel used to display the triangle and receive mouse events.
     */
    private class TrianglePanel extends javax.swing.JPanel {

        public Dimension getPreferredSize() {
            FontMetrics fm = getFontMetrics(getFont());
            int width = sideLength + fm.stringWidth(symbols[0]) + fm.stringWidth(symbols[2]);
            int height = triangle.height + fm.getHeight();
            return new Dimension(width+10, height+10);
        }
            
        public void paintComponent(java.awt.Graphics g) {
            super.paintComponent(g);
            g.drawPolygon(triangle);
            
            FontMetrics fm = getFontMetrics(getFont());
            int textX = triangle.v[0][0] - fm.stringWidth(symbols[0]);
            int textY = triangle.v[0][1];
            g.drawString(symbols[0], textX, textY);
            textX = triangle.v[1][0];
            textY = triangle.v[1][1];
            g.drawString(symbols[1], textX, textY);
            textX = triangle.v[2][0] - fm.stringWidth(symbols[1])/2;
            textY = triangle.v[2][1];//+ fm.getHeight();
            g.drawString(symbols[2], textX, textY);
            
            for(int k=0; k<triangle.nTick; k++) {
                g.drawLine(triangle.ticks0[k][0],triangle.ticks0[k][1],
                            triangle.ticks1[k][0],triangle.ticks1[k][1]);
                g.drawLine(triangle.ticks0[triangle.nTick-k-1][0],triangle.ticks0[triangle.nTick-k-1][1],
                            triangle.ticks2[k][0],triangle.ticks2[k][1]);
                g.drawLine(triangle.ticks1[k][0],triangle.ticks1[k][1],
                            triangle.ticks2[k][0],triangle.ticks2[k][1]);
            }
        }//end paintComponent method
    }//end of TrianglePanel class
    
    //listener for mouse events.  Reads mouse position and converts
    //to fractions
    private class ClickListener extends java.awt.event.MouseAdapter {
        public void mousePressed(java.awt.event.MouseEvent evt) {
            int x = evt.getX();
            int y = evt.getY();
            double[] frac = triangle.fractions(x,y);
            for(SimulationListener.Linker linker=listenerManager.first(); linker!=null; linker=linker.next) {
                ((Listener)linker.listener).ternaryAction(frac[0], frac[1], frac[2]);
            }
            //System.out.println(frac[0]+" "+frac[1]+" "+frac[2]);
        }
    }//end of ClickListener class
    
    /**
     * The height of the triangle, equal to sqrt(3)/2 times the length of a side.
     */
     //placed here rather than in Triangle class so it can be used in Triangle's
     //constructor.  Cannot be declared static in Triangle because Triangle is member class
    private static int h(int L) {return (int)(0.5*Math.sqrt(3.)*L);}

    /**
     * Class used to draw and perform calculations for the triangle.
     */
    private class Triangle extends java.awt.Polygon {
        
        public int height;
        private int x0, y0;
        private final int nTick = 4;
        private double[] x = new double[3];
        public final int[][] v = new int[3][2];//coordinates of vertices v[vertexNo][xy]
        public final int[][] ticks0 = new int[nTick][2];
        public final int[][] ticks1 = new int[nTick][2];
        public final int[][] ticks2 = new int[nTick][2];
        public Triangle(int L, int x0, int y0) {
            super(new int[] {   x0+0,    x0+L, x0+L/2}, //x points       2
                  new int[] {y0+h(L), y0+h(L),   y0+0}, //y points      0 1
                  3);
            height = h(L);
            this.x0 = x0; this.y0 = y0;
            v[0][0] = x0+0;   v[0][1] = y0+height;
            v[1][0] = x0+L;   v[1][1] = y0+height;
            v[2][0] = x0+L/2; v[2][1] = y0+0;
            
            for(int k=0; k<nTick; k++) {//which tick on line
                for(int j=0; j<2; j++) {//x, y
                    ticks0[k][j] = v[0][j] + (k+1)*(v[1][j]-v[0][j])/(nTick+1);
                    ticks1[k][j] = v[0][j] + (k+1)*(v[2][j]-v[0][j])/(nTick+1);
                    ticks2[k][j] = v[1][j] + (k+1)*(v[2][j]-v[1][j])/(nTick+1);
                }
            }
        }
        
        
        /**
         * Returns the ternary fractions corresponding to the given pixel-point 
         * in the triangle.
         */
        public double[] fractions(int xP, int yP) {
            double xx = 0.5*Math.sqrt(3.)*(double)(xP-x0)/(double)height;
            double yy = 0.5*(1.0 - (double)(yP-y0)/(double)height);
            x[0] = 1.0 - xx - yy;
            x[1] = xx - yy;
            x[2] = 2.*yy;
            //set to zero components if outside triangle
            if(x[0] < 0.0) x[0] = 0.0;
            if(x[1] < 0.0) x[1] = 0.0;
            if(x[2] < 0.0) x[2] = 0.0;
            //enforce precision
            x[0] = Math.round((float)(x[0]*power))/(double)power;
            x[1] = Math.round((float)(x[1]*power*10))/(double)(power*10);
            x[2] = Math.round((float)(x[2]*power))/(double)power;
            //renormalize to for case where outside triangle
            double norm = x[0] + x[1] + x[2];
            for(int i=0; i<3; i++) {x[i] /= norm;}
            //update display boxes if present
            if(boxes != null) {
                for(int i=0; i<3; i++) {boxes[i].setText(DisplayBox.format(x[i],precision));};
            }
            return x;
        }//end of fractions
    }//end of Triangle
    
    public interface Listener extends SimulationListener {
        public void ternaryAction(double x1, double x2, double x3);
    }
    
//    public static void main(String[] args) {
//        SimulationGraphic sim = new etomica.simulations.HSMD2D();
//        Simulation.instance = sim;
//        DeviceTernarySelector selector = 
//            new DeviceTernarySelector(sim, new String[] {"Benzene","Toluene","Xylene"});
//        selector.addListener(new DeviceTernarySelector.Listener() {
//            public void ternaryAction(double x1, double x2, double x3) {
//                System.out.println("Event!: " + x1 + " " + x2 + " " + x3);
//            }
//        });
//        selector.setShowValues(true);
//        selector.setSymbols(new String[] {"B", "T", "X"});
//      //  selector.setSideLength(400);
//        sim.elementCoordinator.go();
//        sim.makeAndDisplayFrame();
//    }
    
}