package etomica.virial.cluster2.ui.wizards;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A model represents the domain-specific data on which the application
 * operates. The domain-specific validation of the model is done at the isValid
 * method and operates on a state-by-state basis.
 * 
 * @author Demian Lessa
 */
public class Model {

  private String modelName = null;
  private Map<String, Object> model = null;

  private Model() {

    model = new HashMap<String, Object>();
  }

  private Model(final String name) {

    this();
    modelName = name;
  }

  private void ensureState(String key) {

    if (!model.containsKey(key)) {
      throw new NonexistentStateException(modelName, key);
    }
  }

  private void registerState(final List<String> keys) {

    for (String key : keys) {
      model.put(key, null);
    }
  }

  private void setDefaultModel(Map<String, Object> defaultModel) {

    model.putAll(defaultModel);
  }

  protected boolean isValid(final String key, final Object value) {

    return false;
  }

  public static Model createModel(final String modelName,
      final List<String> keys) {

    Model m = new Model(modelName);
    m.registerState(keys);
    return m;
  }

  public static Model createModel(final String modelName,
      final Map<String, Object> defaultModel) {

    Model m = new Model(modelName);
    m.setDefaultModel(defaultModel);
    return m;
  }

  public final Object getState(final String key) {

    ensureState(key);
    return model.get(key);
  }

  public final void setState(final String key, final Object value) {

    ensureState(key);
    if (!isValid(key, value)) {
      throw new InvalidStateException(modelName, key);
    }
    model.put(key, value);
  }

  class NonexistentStateException extends RuntimeException {

    private static final long serialVersionUID = 8986512369447774563L;

    NonexistentStateException(final String modelName, final String key) {

      super("The key " + key + " does not exist in the model " + modelName);
    }
  }

  class InvalidStateException extends RuntimeException {

    private static final long serialVersionUID = -2294817079959667990L;

    InvalidStateException(final String modelName, final String key) {

      super("The value passed to " + key + " puts the model " + modelName
          + " in an invalid state.");
    }
  }
}