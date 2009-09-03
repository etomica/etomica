package etomica.virial.cluster2.ui.wizards;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import etomica.virial.cluster2.ui.View;
import etomica.virial.cluster2.ui.ViewFactory;

public class Controller {

  public static final int EVENT_INITIALIZE = 0x0000;
  public static final int EVENT_NEXT_VIEW = 0x0001;
  public static final int EVENT_PRIOR_VIEW = 0x0002;
  public static final int EVENT_FINISH = 0x0004;
  public static final int EVENT_CANCEL = 0x0008;
  private List<Model> model = null;
  private ModelFactory modelFactory = null;

  public Controller(final ModelFactory mfactory) {

    this();
    modelFactory = mfactory;
    dispatch(EVENT_INITIALIZE, null);
  }

  private Controller() {

    model = new ArrayList<Model>();
  }

  protected void dispatch(int eventID, final Object param) {

    switch (eventID) {
      case EVENT_INITIALIZE: {
        nextState();
        break;
      }
      case EVENT_NEXT_VIEW: {
        View view = (View) param;
        if (updateCurrentModel(view.getData())) {
          nextState();
        }
        break;
      }
      case EVENT_PRIOR_VIEW: {
        dropCurrentModel();
        nextState();
        break;
      }
      case EVENT_FINISH: {
        View view = (View) param;
        if (updateCurrentModel(view.getData())) {
          finish();
        }
        break;
      }
      case EVENT_CANCEL: {
        cancel();
        break;
      }
      default:
        break;
    }
  }

  public void cancel() {

    // TODO Auto-generated method stub
  }

  public void finish() {

    // TODO Auto-generated method stub
  }

  protected void dropCurrentModel() {

    if (model.size() > 0) {
      model.remove(model.size() - 1);
    }
  }

  protected Model getCurrentModel() {

    if (model.size() == 0) {
      return null;
    }
    return model.get(model.size() - 1);
  }

  protected void nextState() {

    Model nextModel = modelFactory.createNext(getCurrentModel());
    model.add(nextModel);
    View nextView = ViewFactory.createView(getCurrentModel());
    nextView.activate();
  }

  protected boolean updateCurrentModel(Map<String, Object> data) {

    // TODO Auto-generated method stub
    return false;
  }
}