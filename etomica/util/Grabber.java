/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.image.BufferedImage;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;

import etomica.graphics.SimulationGraphic;

public class Grabber {
    public static class CornerData {
        public double x, y;
        public int pixelX = -1, pixelY = -1;
    }
    
    protected final CornerData[] cornerData;
    protected final JFrame f;
    protected final DefaultMouseListener defaultMouseListener;
    protected final ImagePanel imgPanel;
    protected final JButton buttonCorner1, buttonCorner2;
    protected ActionCornerSetter activeCornerSetter;
    protected final JTextArea outputArea;
    
    public static final class ActionCornerSetter implements Action {
        protected MouseListener myListener;
        protected final Grabber grabber;
        protected final int iCorner;
        
        public ActionCornerSetter(Grabber g, int iCorner) {
            grabber = g;
            this.iCorner = iCorner;
        }

        public void actionPerformed(ActionEvent e) {
            if (grabber.activeCornerSetter == this) {
                grabber.activeCornerSetter = null;
                grabber.imgPanel.removeMouseListener(myListener);
                myListener = null;
                grabber.defaultMouseListener.setEnabled(true);
            }
            else {
                grabber.activeCornerSetter = this;
                grabber.defaultMouseListener.setEnabled(false);
                myListener = new CornerSetListener(grabber, this, iCorner);
                grabber.imgPanel.addMouseListener(myListener);
            }
        }
        
        public void deactivate() {
            if (grabber.activeCornerSetter != this) {
                return;
            }
            grabber.activeCornerSetter = null;
            grabber.imgPanel.removeMouseListener(myListener);
            myListener = null;
            grabber.defaultMouseListener.setEnabled(true);
        }            

        public void setEnabled(boolean b) {}

        public void removePropertyChangeListener(PropertyChangeListener listener) {}

        public void putValue(String key, Object value) {}

        public boolean isEnabled() {return true;}

        public Object getValue(String key) {
            return null;
        }

        public void addPropertyChangeListener(PropertyChangeListener listener) {}
    }

    public static final class CornerSetListener implements MouseListener {
        protected final Grabber grabber;
        protected final ActionCornerSetter action;
        protected final int iCorner;
        
        public CornerSetListener(Grabber g, ActionCornerSetter a, int iCorner) {
            grabber = g;
            action = a;
            this.iCorner = iCorner;
        }

        public void mouseReleased(MouseEvent e) {}

        public void mousePressed(MouseEvent e) {}

        public void mouseExited(MouseEvent e) {}

        public void mouseEntered(MouseEvent e) {}

        public void mouseClicked(MouseEvent e) {
            if (e.getButton() != MouseEvent.BUTTON1) {
                return;
            }
            
            action.deactivate();
            
            final int xClick = e.getX();
            final int yClick = e.getY();
            final JDialog dialog = new JDialog(grabber.f, "Enter Corner 1 Values", true);
            
            final JButton okButton = new JButton();

            dialog.setMinimumSize(new Dimension(200, 100));
            JPanel panel = new JPanel(new GridBagLayout());
            panel.setMinimumSize(new Dimension(200, 100));
            GridBagConstraints gbc = new GridBagConstraints();
            gbc.gridx = 0;
            gbc.gridy = 0;
            JLabel xLabel = new JLabel("x value");
            panel.add(xLabel, gbc);
            KeyListener keyListener = new KeyListener() {

                public void keyPressed(KeyEvent ev) {
                }

                public void keyReleased(KeyEvent ev) {
                    if (ev.getKeyCode() == KeyEvent.VK_ENTER) {
                        okButton.doClick();
                    }
                }

                public void keyTyped(KeyEvent ev) {
                }
            };
            final JTextField xField = new JTextField(9);
            xField.addKeyListener(keyListener);
            gbc.gridx = 1;
            panel.add(xField, gbc);
            JLabel yLabel = new JLabel("y value");
            gbc.gridx = 0;
            gbc.gridy = 1;
            panel.add(yLabel, gbc);
            final JTextField yField = new JTextField(9);
            yField.addKeyListener(keyListener);
            gbc.gridx = 1;
            panel.add(yField, gbc);
            
            if (grabber.cornerData[iCorner].pixelX > -1) {
                xField.setText(""+grabber.cornerData[iCorner].x);
                yField.setText(""+grabber.cornerData[iCorner].y);
            }
            
            JPanel buttonPanel = new JPanel(new GridLayout(1, 2));
            gbc.gridx = 0;
            gbc.gridy = 2;
            gbc.gridwidth = 2;
            panel.add(buttonPanel, gbc);
            gbc.gridwidth = 1;

            okButton.setAction(new Action() {
                public void addPropertyChangeListener(
                        PropertyChangeListener listener) {}

                public Object getValue(String key) {return null;}

                public boolean isEnabled() {return true;}

                public void putValue(String key, Object value) {}

                public void removePropertyChangeListener(
                        PropertyChangeListener listener) {}

                public void setEnabled(boolean b) {}

                public void actionPerformed(ActionEvent ev) {
                    double xValue, yValue;
                    try {
                        xValue = Double.parseDouble(xField.getText());
                        yValue = Double.parseDouble(yField.getText());
                    }
                    catch (NumberFormatException ex) {
                        // do nothing
                        System.out.println("nope");
                        return;
                    }
                    if (Double.isInfinite(xValue) || Double.isNaN(xValue) || Double.isInfinite(yValue) || Double.isNaN(yValue)) {
                        System.out.println("nope #2");
                        return;
                    }
                    CornerData cd = grabber.cornerData[iCorner];
                    cd.x = xValue;
                    cd.y = yValue;
                    cd.pixelX = xClick;
                    cd.pixelY = yClick;
                    dialog.setVisible(false);
                    dialog.dispose();
                    
                    grabber.outputArea.setText("");
                }
            });
            okButton.setText("OK");
            buttonPanel.add(okButton);

            JButton cancelButton = new JButton();
            cancelButton.setAction(new Action() {
                public void addPropertyChangeListener(
                        PropertyChangeListener listener) {}

                public Object getValue(String key) {return null;}

                public boolean isEnabled() {return true;}

                public void putValue(String key, Object value) {}

                public void removePropertyChangeListener(
                        PropertyChangeListener listener) {}

                public void setEnabled(boolean b) {}

                public void actionPerformed(ActionEvent ev) {
                    dialog.setVisible(false);
                    dialog.dispose();
                }
            });
            cancelButton.setText("Cancel");
            gbc.gridx = 1;
            gbc.gridy = 2;
            buttonPanel.add(cancelButton);
            
            dialog.add(panel);
            dialog.setVisible(true);
        }
    }

    public static final class DefaultMouseListener implements MouseListener {
        protected final Grabber grabber;
        
        public DefaultMouseListener(Grabber g) {
            grabber = g;
        }
        
        public void mouseReleased(MouseEvent e) {}
        public void mousePressed(MouseEvent e) {}
        public void mouseExited(MouseEvent e) {}
        public void mouseEntered(MouseEvent e) {}

        public void mouseClicked(MouseEvent e) {
            if (!isEnabled || e.getButton() != MouseEvent.BUTTON1) {
                return;
            }
            int xPixelClick = e.getX();
            int yPixelClick = e.getY();
            if (grabber.cornerData[0].pixelX < 0 || grabber.cornerData[1].pixelX < 1) {
                return;
            }
            CornerData c0 = grabber.cornerData[0];
            double x0 = c0.x;
            double y0 = c0.y;
            int px0 = c0.pixelX;
            int py0 = c0.pixelY;
            CornerData c1 = grabber.cornerData[1];
            double x1 = c1.x;
            double y1 = c1.y;
            int px1 = c1.pixelX;
            int py1 = c1.pixelY;
            
            double xClick = x0 + (x1-x0)*(xPixelClick-px0)/(px1-px0);
            double yClick = y0 + (y1-y0)*(yPixelClick-py0)/(py1-py0);
            grabber.outputArea.append(xClick+" "+yClick+"\n");
//            System.out.println(xClick+" "+yClick);
        }
        
        public void setEnabled(boolean newIsEnabled) {
            isEnabled = newIsEnabled;
        }
        
        public boolean getEnabled() {return isEnabled;}
        
        protected boolean isEnabled = true;
    }

    public static void main(String[] args) {
        String filename = "/tmp/plot.png";
        if (args.length > 0) {
            filename = args[0];
        }
        new Grabber(filename);
    }
     
    public Grabber(String filename) {
        cornerData = new CornerData[2];
        cornerData[0] = new CornerData();
        cornerData[1] = new CornerData();
        
        BufferedImage img;
        try {
            img = ImageIO.read(new File(filename)); //args[0]));
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        int imgHeight = img.getHeight();
        int imgWidth = img.getWidth();

        f = new JFrame();
        f.setSize(imgWidth+ 50,imgHeight+300);
        f.pack();
        f.setVisible(true);
        f.addWindowListener(SimulationGraphic.WINDOW_CLOSER);
        
        JPanel topPanel = new JPanel(new GridBagLayout());
        f.getContentPane().add(topPanel);
        
        JPanel imgContainer = new JPanel(new GridLayout(1,1));
        
        //cp.addMouseListener(new MouseEventClusterPanel());
        Dimension d = new Dimension();
        d.width=imgWidth;
        d.height=imgHeight;
        imgContainer.setPreferredSize(d);
        imgContainer.setMinimumSize(d);
        imgContainer.setMaximumSize(d);
        imgContainer.setBackground(null);
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        topPanel.add(imgContainer, gbc);

        
        imgPanel = new ImagePanel();
        defaultMouseListener = new DefaultMouseListener(this);
        imgPanel.addMouseListener(defaultMouseListener);
        imgContainer.add(imgPanel);
        imgPanel.setImage(img);
        
        JPanel controlPanel = new JPanel(new GridBagLayout());
        JPanel cornerButtonPanel = new JPanel(new GridLayout(3, 1));
        buttonCorner1 = new JButton(new ActionCornerSetter(this, 0));
        buttonCorner1.setText("Set Corner #1");
        cornerButtonPanel.add(buttonCorner1);
        buttonCorner2 = new JButton(new ActionCornerSetter(this, 1));
        buttonCorner2.setText("Set Corner #2");
        cornerButtonPanel.add(buttonCorner2);
        JButton clearButton = new JButton(new Action() {

            public void addPropertyChangeListener(
                    PropertyChangeListener listener) {}

            public Object getValue(String key) {return null;}

            public boolean isEnabled() {return true;}

            public void putValue(String key, Object value) {}

            public void removePropertyChangeListener(
                    PropertyChangeListener listener) {}

            public void setEnabled(boolean b) {}

            public void actionPerformed(ActionEvent e) {
                outputArea.setText("");
            }
        });
        clearButton.setText("Clear");
        cornerButtonPanel.add(clearButton);
        
        gbc.gridx = 0;
        gbc.gridy = 0;
        controlPanel.add(cornerButtonPanel, gbc);
        
        outputArea = new JTextArea(10, 40);
        JScrollPane outputScrollPane = new JScrollPane(outputArea);
        outputArea.setEditable(false);
        gbc.gridx = 1;
        controlPanel.add(outputScrollPane, gbc);
        
        gbc.gridx = 0;
        gbc.gridy = 1;
        topPanel.add(controlPanel, gbc);
    }

    public static class ImagePanel extends Component {

        /**
         * Constructor sets default size of the plot
         */
        public ImagePanel () {
            this(Color.white,400,400);
        }

        public ImagePanel(Color background, int xSize, int ySize) {
            setBackground(background);
            setSize(xSize,ySize);
        }

        public void setImage(Image newImage) {
            myImage = newImage;
        }

        public void update(Graphics g) {
            paint(g);
        }

        public void paint (Graphics g) {
            g.drawImage(myImage, 0, 0, null);
        }

        protected Image myImage;
        protected Graphics osg;
    }
}
