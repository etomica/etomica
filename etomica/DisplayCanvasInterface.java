package etomica;
import java.awt.*;
import java.awt.event.*;

public interface DisplayCanvasInterface {
    public void createOffScreen();
    public void createOffScreen(int p);
    public void createOffScreen(int w, int h);
    public void doPaint(Graphics g);
    public void update(Graphics g);
    public void paint(Graphics g);
    public void setPhase(Phase p);
    public void setMovable(boolean b);
    public boolean isMovable();
    public void setResizable(boolean b);
    public boolean isResizable();
    public void setWriteScale(boolean s);
    public boolean getWriteScale();
    public void setHighQuality(boolean q);
    public boolean getHighQuality();
    public void initialize();
    public void addMouseListener(MouseListener listener);
    public void addMouseMotionListener(MouseMotionListener listener);
    public void addKeyListener(KeyListener listener);
    public void setSize(int width, int height);
    public void setMinimumSize(Dimension temp);
    public void setMaximumSize(Dimension temp);
    public void setPreferredSize(Dimension temp);
    public void setBounds(int x, int y, int width, int height);
    public java.awt.Dimension getSize();
    public void repaint();
    public void requestFocus();
    public void transferFocus();
    public void setShiftX(float x);
    public void setShiftY(float y);
    public void setPrevX(float x);
    public void setPrevY(float y);
    public void setXRot(float x);
    public void setYRot(float y);
    public void setZoom(float z);
    public float getShiftX();
    public float getShiftY();
    public float getPrevX();
    public float getPrevY();
    public float getXRot();
    public float getYRot();
    public float getZoom();
}

