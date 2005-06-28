package etomica.graphics;
import java.awt.Color;
import java.awt.Graphics;

import etomica.data.DataSourceScalar;
import etomica.data.types.DataDouble;

    public class DisplayScrollingGraph extends Display {

        Color bgColor = new Color(255, 255, 240);
        int margin = 10;
        int number = 0;
        int current = 0;
        int pixels = 300;
        int stepsToGraph = 100;
        boolean useMeterStatistics;
        double[] values = new double[stepsToGraph];
        double vMin = 0.0;
        double vMax = 60.0;
        String legend = new String("");
        int ticks = 5;
        double[] tickValues = {0.0, 15.0, 30.0, 45.0, 60.0};
        double xScale = (pixels - 2 * margin) / (double) stepsToGraph;
        double yScale = (pixels - 2 * margin) / (vMax - vMin);
        double vAvg;
        double v2Avg;
        double stdDev;
        DataSourceScalar meter;
        etomica.units.Unit unit;
        
        public DisplayScrollingGraph() {
            useMeterStatistics = false;
            clear();
        }
        
        void clear () {
            current = number = 0;
            for (int v = 1; v < values.length; v++)
                values[v] = 0;
        }

        void addValue(double value) {
            if (values.length < stepsToGraph) {
                double[] temp = values;
                values = new double[stepsToGraph];
                for (int v = 0; v < temp.length; v++)
                    values[v] = temp[v];
            }
            if (number < stepsToGraph - 1)
                ++number;
            if (++current == stepsToGraph)
                current = 0;
            values[current] = value;
        }
        
    /**
     * Accessor method for the meter that generates the displayed value.
     */
    public DataSourceScalar getMeter() {
        System.out.println("DisplayBox.getMeter value = "+meter);
        return meter;
    }

    /**
     * Specifies the meter that generates the displayed value.
     */
    public void setMeter(DataSourceScalar m) {
        meter = m;
        setUnit(m.getDataInfo().getDimension().defaultIOUnit());
    }
    
    
    /**
     * Accessor method to set the physical units of the displayed value.
     * Text describing unit is used in label.
     */
    public void setUnit(etomica.units.Unit u) {
        unit = u;
//        setLabel();
    }

    /**
     * Returns the physical units of the displayed value.
     */
    public etomica.units.Unit getUnit() {return unit;}

    public void doUpdate() {
        addValue(((DataDouble)meter.getData()).x); //re-examine this
    }
        
    private void updateStatistics() {
        if(useMeterStatistics) {
            vAvg = meter.average();
            stdDev = Math.sqrt(meter.variance());
        }
        else {
            vAvg = v2Avg = 0;
            for (int n = 0; n < number; n++) {
                vAvg += values[n];
                v2Avg += values[n] * values[n];
            }
            stdDev = 0;
            if (number > 0) {
                vAvg /= number;
                v2Avg /= number;
                stdDev = Math.sqrt(v2Avg - vAvg * vAvg);
            }
        }
    }
    
    public void resetScale() {          
        double min = +Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;
        for(int i=0; i<values.length; i++) {
            min = (values[i] < min) ? values[i] : min;
            max = (values[i] > max) ? values[i] : max;
        }
        setVMin(min-1.1*stdDev);
        setVMax(max+1.1*stdDev);
    }
    
    public void paint(Graphics g) {doPaint(g);}
    public void doPaint(Graphics osg) {

        osg.setColor(bgColor);
        osg.fillRect(0, 0, pixels, pixels);
        updateStatistics();
        int x = margin + (int) (number * xScale);
        int y = margin + (int) ((vMax - vAvg) * yScale);
        int dy = (int) (stdDev * yScale);
        osg.setColor(Color.lightGray);
        osg.fillRect(margin, y - dy, x - margin, 2 * dy);
        osg.setColor(Color.magenta);
        osg.drawLine(margin, y, x, y);

        osg.setColor(Color.blue);
        int xOld = 0;
        int yOld = 0;
        for (int n = 0; n < number; n++) {
            x = margin + (int) (n * xScale);
            y = current - number + n;
            if (y < 0)
                y += stepsToGraph;
            y = margin + (int) ((vMax - values[y]) * yScale);
            if (n == 0)
                osg.drawLine(x, y, x, y);
            else osg.drawLine(xOld, yOld, x, y);
            xOld = x;
            yOld = y;
        }

        x = margin;
        if (tickValues != null) {
            for (int tick = 0; tick < ticks; tick++) {
                y = margin + (int) ((vMax - tickValues[tick]) * yScale);
                osg.setColor(Color.cyan);
                osg.drawLine(margin, y, pixels - margin, y);
                y += margin / 2;
                osg.setColor(Color.black);
                osg.drawString("" + tickValues[tick], x, y);
            }
        }
        osg.setColor(Color.red);
        osg.drawString(legend, pixels / 4, pixels - margin / 2);
    }

    public void setNTicks(int n) {
        ticks = n;
        resetTickValues();
    }
    public int getNTicks() {return ticks;}
    
    public void setVMin(double v) {
        vMin = v;
        yScale = (pixels - 2 * margin) / (vMax - vMin);
        resetTickValues();
    }
    public double getVMin() {return vMin;}
    public void setVMax(double v) {
        vMax = v;
        yScale = (pixels - 2 * margin) / (vMax - vMin);
        resetTickValues();
    }
    public double getVMax() {return vMax;}
    
    private void resetTickValues() {
        double[] newTicks = new double[ticks];
        double delta = (vMax - vMin)/(ticks-1);
        for(int i=0; i<ticks; i++) newTicks[i] = vMin + i*delta;
        tickValues = newTicks;
    }
            
    
    public void setUseMeterStatistics(boolean b) {useMeterStatistics = b;}
    public boolean getUseMeterStatistics() {return useMeterStatistics;}
}

