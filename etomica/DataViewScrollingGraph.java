package simulate;
import java.awt.*;

    public class DataViewScrollingGraph extends DataView {

        Color bgColor = new Color(255, 255, 240);
        int margin = 10;
        int number = 0;
        int current = 0;
        int stepsToGraph = 100;
        double[] values = new double[stepsToGraph];
        double vMin = -1000;
        double vMax = 1000;
        String legend = new String("");
        int ticks = 5;
        double[] tickValues = {-1000.0, -500, 0.0, 500, 1000};

        void clear () {
            current = number = 0;
            for (int v = 1; v < values.length; v++)
                values[v] = 0;
        }

        void add (double value) {
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

        public void paint (Graphics g) {
            update(g);
        }

        public void update (Graphics osg) {
            add(parentDisplay.firstMeter.currentValue(parentDisplay.phase));
            osg.setColor(bgColor);
            int pixels = parentDisplay.pixels;
            osg.fillRect(0, 0, pixels, pixels);

            double vAvg = 0;
            double v2Avg = 0;
            for (int n = 0; n < number; n++) {
                vAvg += values[n];
                v2Avg += values[n] * values[n];
            }
            double stdDev = 0;
            if (number > 0) {
                vAvg /= number;
                v2Avg /= number;
                stdDev = Math.sqrt(v2Avg - vAvg * vAvg);
            }
            double xScale = (pixels - 2 * margin) / (double) stepsToGraph;
            double yScale = (pixels - 2 * margin) / (vMax - vMin);
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

    public void setVMin(double v) {vMin = v;}
    public void setVMax(double v) {vMax = v;}
    }

