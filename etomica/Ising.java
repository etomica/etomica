// ComPhys File: Ising.java
// Chapter 17: Monte Carlo simulation of the Ising model on
//             a square lattice using the Metropolis Algorithm

import java.applet.*;
import java.awt.*;
import java.awt.event.*;
import comphys.*;

public class Ising extends Applet implements ActionListener, ItemListener, Runnable {

    // Ising model quantities

    int L = 16;                 // linear dimension of lattice
    int N = L * L;              // number of spins
    int[][] spin;               // 2-d array of spin variables
    double J = 1;               // exchange coupling
    double T = 2;               // temperature
    double H = 0;               // magnetic field
    double[][] w;               // Boltzmann factors

    double E;                   // total energy
    int M;                      // total magnetic moment
    int SM;                     // staggered magnetization
    int SS;                     // sum s_i * s_j
    boolean hotStart = true;    // random initial configuration

    void initial () {

        // allocate spin array if necessary
        if (spin == null || spin.length < L)
            spin = new int[L][L];
        N = L * L;

        // initialize spins
        M = 0;
        SM = 0;
        for (int x = 0; x < L; x++) {
            for (int y = 0; y < L; y++) {

                spin[x][y] = 1;
                if (hotStart)
                    if (Math.random() > 0.5)
                        spin[x][y] = -1;
                M += spin[x][y];
                if ( (x + y) % 2 == 0 )
                    SM += spin[x][y];
                else
                    SM -= spin[x][y];

            }
        }

        // compute total energy
        SS = 0;
        for (int x = 0; x < L; x++) {

            int right = x + 1;
            if (right == L)
                right = 0;
            for (int y = 0; y < L; y++) {

                int up = y + 1;
                if (up == L)
                    up = 0;
                int sum = spin[x][up] + spin[right][y];
                SS += spin[x][y] * sum;

            }
        }

        E = - J * SS - H * M;

        computeBoltzmannFactors();
        initializeSums();

    }

    void computeBoltzmannFactors () {

        // compute Boltzmann probability ratios
        if (w == null)
            w = new double[17][3];

        for (int i = -8; i <= 8; i += 4) {

            w[i + 8][0] = Math.exp( - (i * J + 2 * H) / T);
            w[i + 8][2] = Math.exp( - (i * J - 2 * H) / T);

        }

    }

    int sumOfNeighbors (int x, int y) {

        // periodic boundary conditions

        int left = 0;
        if (x == 0)
            left = spin[L - 1][y];
        else left = spin[x - 1][y];
        int right = 0;
        if (x == L - 1)
            right = spin[0][y];
        else right = spin[x + 1][y];
        int down = 0;
        if (y == 0)
            down = spin[x][L - 1];
        else down = spin[x][y - 1];
        int up = 0;
        if (y == L - 1)
            up = spin[x][0];
        else up = spin[x][y + 1];

        return left + right + up + down;

    }

    int mcs;                    // Monte Carlo steps per spin
    int accept;                 // accepted Metropolis steps

    void Metropolis () {

        // one Monte Carlo step per spin
        for (int ispin = 0; ispin < N; ispin++) {

            int x = (int) (L * Math.random());
            if (x == L)
                --x;
            int y = (int) (L * Math.random());
            if (y == L)
                --y;
            int dSS = 2 * spin[x][y] * sumOfNeighbors(x, y);
            if (Math.random() < w[dSS + 8][spin[x][y] + 1]) {

                spin[x][y] = -spin[x][y]; // flip spin
                ++accept;
                M += 2 * spin[x][y];
                if ( (x + y) % 2 == 0 )
                    SM += 2 * spin[x][y];
                else
                    SM -= 2 * spin[x][y];
                SS -= dSS;

            }

        }

        ++mcs;
        data();

    }

    double Esum;                // accumulator to compute <E>
    double EEsum;               // accumulator to compute <E*E>
    double Msum;                // accumulator to compute <M>
    double MMsum;               // accumulator to compute <M*M>
    double absMsum;             // accumulator to compute <|M|>

    void initializeSums () {

        Esum = 0;
        EEsum = 0;
        Msum = 0;
        MMsum = 0;
        absMsum = 0;
        mcs = 0;
        accept = 0;

    }

    double ePerSpin;
    double mPerSpin;

    void data () {

        E = - J * SS - H * M;
        ePerSpin = E / N;
        Esum += E;
        EEsum += E * E;

        if (J < 0) {            // anti-ferromagnetic case

            mPerSpin = SM / (double) N;
            Msum += SM;
            MMsum += SM * (double)SM;
            absMsum += Math.abs(SM);

        } else {

            mPerSpin = M / (double) N;
            Msum += M;
            MMsum += M * (double)M;
            absMsum += Math.abs(M);

        }

    }

    double acceptAverage;
    double eAverage;
    double e2Average;
    double mAverage;
    double m2Average;
    double abs_mAverage;
    double cPerSpin;
    double chiPerSpin;

    void computeAverages () {

        double norm = 1.0 / N;
        if (mcs > 0)
            norm /= mcs;
        acceptAverage = accept * norm;
        eAverage = Esum * norm;
        e2Average = EEsum * norm;
        mAverage = Msum * norm;
        m2Average = MMsum * norm;
        abs_mAverage = absMsum * norm;
        cPerSpin = (e2Average - N * eAverage * eAverage) / (T * T);
        chiPerSpin = (m2Average - N * mAverage * mAverage) / T;

    }

    class Picture extends Canvas {

        int pixels = 200;

        Picture () {

            setSize(pixels, pixels);
            setBackground(Color.white);

        }

        Image offScreen;
        Graphics osg;

        void createOffScreen () {

            if (offScreen == null) {

                offScreen = createImage(pixels, pixels);
                osg = offScreen.getGraphics();

            }

        }

    }

    class Lattice extends Picture {

        public void paint (Graphics g) {

            update(g);

        }

        public void update (Graphics g) {

            createOffScreen();
            osg.setColor(Color.black);
            osg.fillRect(0, 0, pixels, pixels);

            int d = pixels / L;
            int margin = (pixels - d * L) / 2;
            int sep = 1;
            if (d < 3)
                sep = 0;

            for (int x = 0; x < L; x++) {
                for (int y = 0; y < L; y++) {

                    if (spin[x][y] == 1)
                        osg.setColor(Color.red);
                    else osg.setColor(Color.green);
                    osg.fillRect(x * d + margin, y * d + margin,
                                 d - sep, d - sep);

                }
            }

            g.drawImage(offScreen, 0, 0, null);

        }

    }

    int stepsToGraph = 100;

    class ScrollingGraph extends Picture {

        Color bgColor = new Color(255, 255, 240);
        int margin = 10;
        int number = 0;
        int current = 0;
        double[] values = new double[stepsToGraph];
        double vMin = -1;
        double vMax = 1;
        String legend = new String("");
        int ticks;
        double[] tickValues;

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

        public void update (Graphics g) {

            createOffScreen();
            osg.setColor(bgColor);
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

            g.drawImage(offScreen, 0, 0, null);

        }

    }

    class Output extends Picture {

        public void paint (Graphics g) {

            update(g);

        }

        public void update (Graphics g) {

            createOffScreen();
            osg.setColor(Color.yellow);
            osg.fillRect(0, 0, pixels, pixels);

            int x = 10;
            int dx = pixels / 2;
            int y = 0;
            int dy = 20;

            osg.setColor(Color.black);
            osg.drawString("MC steps", x, y += dy);
            osg.drawString("" + mcs, x + dx, y);
            osg.drawString("<acceptance>", x, y += dy);
            osg.drawString(CP.format(acceptAverage, 4), x + dx, y);
            osg.drawString("<E / spin>", x, y += dy);
            osg.drawString(CP.format(eAverage, 4), x + dx, y);
            osg.drawString("<M / spin>", x, y += dy);
            osg.drawString(CP.format(mAverage, 4), x + dx, y);
            osg.drawString("<|M| / spin>", x, y += dy);
            osg.drawString(CP.format(abs_mAverage, 4), x + dx, y);
            osg.drawString("C / spin", x, y += dy);
            osg.drawString(CP.format(cPerSpin, 4), x + dx, y);
            osg.drawString("chi / spin", x, y += dy);
            osg.drawString(CP.format(chiPerSpin, 4), x + dx, y);
            if (J < 0)
                osg.drawString("J < 0 staggered M and chi", x, y += dy);

            g.drawImage(offScreen, 0, 0, null);

        }

    }

    Lattice lattice;
    ScrollingGraph mGraph, eGraph;
    Output output;

    TextField LTextField;
    TextField JTextField;
    TextField HTextField;
    TextField TTextField;
    Checkbox hotBox, coldBox;
    Button resetButton;
    TextField skipTextField;
    int stepsToSkip;
    TextField graphTextField;
    Button thermButton;
    Button runButton;


    public void init () {

        initial();

        Panel picturePanel = new Panel();
        picturePanel.setLayout(new GridLayout(2, 2));
        picturePanel.add(lattice = new Lattice());
        picturePanel.add(mGraph = new ScrollingGraph());
        mGraph.legend = new String("magnetization per spin");
        mGraph.ticks = 5;
        double[] mTicks = {-1.0, -0.5, 0.0, 0.5, 1.0};
        mGraph.tickValues = mTicks;
        picturePanel.add(eGraph = new ScrollingGraph());
        eGraph.legend = new String("energy per spin");
        eGraph.vMin = -2;
        eGraph.vMax = 0;
        eGraph.ticks = 5;
        double[] eTicks = {-2.0, -1.5, -1.0, -0.5, 0.0};
        eGraph.tickValues = eTicks;
        picturePanel.add(output = new Output());

        add(picturePanel);

        Panel p = new Panel();
        p.setLayout(new GridLayout(8, 2));

        p.add(new Label("Lattice size L ="));
        p.add(LTextField = new TextField("" + L));
        LTextField.addActionListener(this);

        p.add(new Label("Coupling J ="));
        p.add(JTextField = new TextField("" + J));
        JTextField.addActionListener(this);

        p.add(new Label("Magnetic Field H ="));
        p.add(HTextField = new TextField("" + H));
        HTextField.addActionListener(this);

        p.add(new Label("Temperature T ="));
        p.add(TTextField = new TextField("" + T));
        TTextField.addActionListener(this);

        Panel pp = new Panel();
        CheckboxGroup cg = new CheckboxGroup();
        pp.add(hotBox = new Checkbox("Hot", cg, hotStart));
        hotBox.addItemListener(this);
        pp.add(coldBox = new Checkbox("Cold", cg, !hotStart));
        coldBox.addItemListener(this);
        p.add(pp);

        p.add(resetButton = new Button("Reset"));
        resetButton.addActionListener(this);

        p.add(new Label("MC steps to skip:"));
        p.add(skipTextField = new TextField("" + stepsToSkip));
        skipTextField.addActionListener(this);

        p.add(new Label("MC steps to graph:"));
        p.add(graphTextField = new TextField("" + stepsToGraph));
        graphTextField.addActionListener(this);

        p.add(runButton = new Button("Start"));
        runButton.addActionListener(this);

        p.add(thermButton = new Button("Junk therm steps"));
        thermButton.addActionListener(this);

        add(p);

    }

    Thread runThread;
    boolean running;

    public void actionPerformed (ActionEvent event) {

        boolean needToReset = false;

        if (event.getSource() == LTextField) {

            int newL = L;

            try {

                newL = Integer.parseInt(event.getActionCommand());

            } catch (NumberFormatException nfe) { }

            if (L != newL && newL > 1) {

                running = false;
                needToReset = true;
                L = newL;

            }

            LTextField.setText("" + L);

        } else if (event.getSource() == JTextField) {

            try {

                J = new Double(event.getActionCommand()).doubleValue();

            } catch (NumberFormatException nfe) { }

            JTextField.setText("" + J);
            computeBoltzmannFactors();
            eGraph.vMin = Math.min(-2, -2 * Math.abs(J) - Math.abs(H));

        } else if (event.getSource() == HTextField) {

            try {

                H = new Double(event.getActionCommand()).doubleValue();

            } catch (NumberFormatException nfe) { }

            HTextField.setText("" + H);
            computeBoltzmannFactors();
            eGraph.vMin = Math.min(-2, -2 * Math.abs(J) - Math.abs(H));

        } else if (event.getSource() == TTextField) {

            try {

                T = new Double(event.getActionCommand()).doubleValue();

            } catch (NumberFormatException nfe) { }

            TTextField.setText("" + T);
            computeBoltzmannFactors();

        } else if (event.getSource() == resetButton) {

            running = false;
            needToReset = true;

        } else if (event.getSource() == skipTextField) {

            try {

                stepsToSkip = Integer.parseInt(event.getActionCommand());

            } catch (NumberFormatException nfe) { }

            skipTextField.setText("" + stepsToSkip);

        } else if (event.getSource() == graphTextField) {

            int newStepsToGraph = stepsToGraph;

            try {

                newStepsToGraph = Integer.parseInt(event.getActionCommand());

            } catch (NumberFormatException nfe) { }

            if (newStepsToGraph != stepsToGraph && newStepsToGraph > 1) {

                stepsToGraph = newStepsToGraph;
                eGraph.clear();
                mGraph.clear();

            }

            graphTextField.setText("" + stepsToGraph);

        } else if (event.getSource() == runButton) {

            if (running) {

                running = false;
                runButton.setLabel("Start");

            } else {

                running = true;
                runThread = new Thread(this);
                runThread.start();
                runButton.setLabel("Stop");

            }

        } else if (event.getSource() == thermButton) {

            initializeSums();
            eGraph.clear();
            mGraph.clear();
            computeAverages();

        }

        if (needToReset) {

            initial();
            eGraph.clear();
            mGraph.clear();
            computeAverages();

        }

        if (!running)
            runButton.setLabel("Start");

        lattice.repaint();
        eGraph.repaint();
        mGraph.repaint();
        output.repaint();

    }

    public void itemStateChanged (ItemEvent itemEvent) {

        hotStart = hotBox.getState();

    }

    void updateGraphData () {

        eGraph.add(ePerSpin);
        mGraph.add(mPerSpin);
        computeAverages();

    }

    public void run () {

        runThread.setPriority(Thread.MIN_PRIORITY);

        while (running) {

            for (int step = 0; step < stepsToSkip; step++) {
                Metropolis();
                updateGraphData();
            }

            Metropolis();
            computeAverages();
            lattice.repaint();
            updateGraphData();
            eGraph.repaint();
            mGraph.repaint();
            output.repaint();

            try {
                runThread.sleep(10);
            } catch (InterruptedException ie) { }

        }

        runThread = null;

    }

    public static void main (String[] args) {

        Ising ising = new Ising();
        CPFrame aFrame = new CPFrame("The Ising Model");

        aFrame.add(ising);
        ising.init();
        aFrame.setSize(700, 450);
        aFrame.setLocation(50, 50);
        aFrame.setVisible(true);

    }

}
