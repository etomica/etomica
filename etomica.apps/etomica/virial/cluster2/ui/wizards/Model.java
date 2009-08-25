package etomica.virial.cluster2.ui.wizards;

import java.util.HashMap;
import java.util.Map;

public class Model {

  private static final Integer MODEL_ID = 0;
  private int modelID = 0;
  private Map<Integer, Object> model = new HashMap<Integer, Object>();

  public Model() {

    model.put(MODEL_ID, ++modelID);
  }

  public final Object checkState(int key) {

    return model.get(key);
  }

  public final boolean modifyState(int key, final Object value) {

    if (!validateState(key, value)) {
      return false;
    }
    model.put(key, value);
    return true;
  }

  protected boolean validateState(int key, final Object value) {

    return (key != MODEL_ID);
  }
}
