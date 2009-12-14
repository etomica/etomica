package etomica.virial.cluster2.mvc.view;

import java.awt.Color;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.KeyAdapter;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.virial.cluster2.mvc.DefaultState;
import etomica.virial.cluster2.mvc.MVCException;
import etomica.virial.cluster2.mvc.State;
import etomica.virial.cluster2.mvc.View;
import etomica.virial.cluster2.mvc.ViewResponse;
import etomica.virial.cluster2.mvc.ViewResponseListener;
import etomica.virial.cluster2.mvc.ViewStatus;
import etomica.virial.cluster2.mvc.WizardPageView;

public class DefaultWizardPage implements WizardPageView {

  private static Boolean onOff = true;
  private ViewResponseListener viewResponseListener;
  private State data = new DefaultState();
  private MouseAdapter mouseAdapter;
  private KeyAdapter keyAdapter;

  public DefaultWizardPage() {

    keyAdapter = new KeyAdapter() {

      @Override
      public void keyPressed(KeyEvent e) {

        if (e.getKeyCode() == KeyEvent.VK_SPACE || e.getKeyCode() == KeyEvent.VK_ENTER) {
          if (e.getComponent() == data.getProperty(ClusterWizard.KEY_NEXT_BUTTON)) {
            doNext();
          }
          else if (e.getComponent() == data.getProperty(ClusterWizard.KEY_BACK_BUTTON)) {
            doBack();
          }
          else if (e.getComponent() == data.getProperty(ClusterWizard.KEY_CANCEL_BUTTON)) {
            doCancel();
          }
          else if (e.getComponent() == data.getProperty(ClusterWizard.KEY_HELP_BUTTON)) {
            doHelp();
          }
          else if (e.getComponent() == data.getProperty(ClusterWizard.KEY_FINISH_BUTTON)) {
            doFinish();
          }
        }
      }
    };

    mouseAdapter = new MouseAdapter() {

      @Override
      public void mouseClicked(MouseEvent e) {

        if (e.getComponent() == data.getProperty(ClusterWizard.KEY_NEXT_BUTTON)) {
          doNext();
        }
        else if (e.getComponent() == data.getProperty(ClusterWizard.KEY_BACK_BUTTON)) {
          doBack();
        }
        else if (e.getComponent() == data.getProperty(ClusterWizard.KEY_CANCEL_BUTTON)) {
          doCancel();
        }
        else if (e.getComponent() == data.getProperty(ClusterWizard.KEY_HELP_BUTTON)) {
          doHelp();
        }
        else if (e.getComponent() == data.getProperty(ClusterWizard.KEY_FINISH_BUTTON)) {
          doFinish();
        }
      }
    };
  }

  public void attach(String key, Object object) {

    data.setProperty(key, object);

    if (object instanceof JButton) {
      ((JButton) object).addMouseListener(mouseAdapter);
      ((JButton) object).addKeyListener(keyAdapter);
    }
    else {
      JPanel panel = (JPanel) object;
      if (panel == data.getProperty(ClusterWizard.KEY_MODEL_PANE)) {
        panel.setBackground(Color.CYAN);
        panel.setOpaque(onOff);
        onOff = !onOff;
        Double randomSeed = Math.random();
        JLabel label = new JLabel("Called from attach: @" + randomSeed);
        panel.add(label);
        if (randomSeed > 0.5) {
          JButton button = new JButton("Click me, dude");
          panel.add(button);
        }
      }
    }
  }

  protected void doFinish() {

    final WizardPageView thisView = this;
    viewResponseListener.onViewResponse(new ViewResponse() {

      public ViewStatus getStatus() {

        return ViewStatus.COMPLETE_SUCCESS;
      }

      public View getView() {

        return thisView;
      }

      public State getData() {

        // TODO Auto-generated method stub
        return null;
      }

      public List<MVCException> getErrors() {

        // TODO Auto-generated method stub
        return null;
      }

    });
  }

  protected void doHelp() {

    // viewResponseListener.onViewResponse(viewResponse);
  }

  protected void doCancel() {

    // viewResponseListener.onViewResponse(viewResponse);
  }

  protected void doBack() {

    // viewResponseListener.onViewResponse(viewResponse);
  }

  protected void doNext() {

    final WizardPageView thisView = this;
    viewResponseListener.onViewResponse(new ViewResponse() {

      public ViewStatus getStatus() {

        return ViewStatus.CONTINUE_NEXT;
      }

      public View getView() {

        return thisView;
      }

      public State getData() {

        // TODO Auto-generated method stub
        return null;
      }

      public List<MVCException> getErrors() {

        // TODO Auto-generated method stub
        return null;
      }

    });
  }

  public void detach(String key, Object object) {

    if (object instanceof JButton) {
      JButton button = (JButton) object;
      button.removeMouseListener(mouseAdapter);
      button.removeKeyListener(keyAdapter);
    }
    if (key.equals(ClusterWizard.KEY_MODEL_PANE)) {
      ((JPanel) object).removeAll();
    }
  }

  public void configure(State state) {

    for (String key : state.getKeys()) {
      data.setProperty(key, state.getProperty(key));
    }
  }

  public void display() {

    JPanel modelPane = (JPanel) data.getProperty(ClusterWizard.KEY_MODEL_PANE);
    modelPane.validate();
    modelPane.repaint();
    JPanel figurePane = (JPanel) data.getProperty(ClusterWizard.KEY_FIGURE_PANE);
    figurePane.validate();
    figurePane.repaint();
  }

  public void setResponseListener(ViewResponseListener listener) {

    this.viewResponseListener = listener;
  }
}