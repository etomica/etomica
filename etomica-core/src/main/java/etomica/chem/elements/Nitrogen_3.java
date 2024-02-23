package etomica.chem.elements;

    public class Nitrogen_3 extends ElementChemical {

        protected Nitrogen_3(String symbol) {
            this(symbol, 14.01);
        }

        protected Nitrogen_3(String symbol, double mass) {
            super(symbol, mass, 7);
        }

        public static final etomica.chem.elements.Nitrogen_3 INSTANCE = new etomica.chem.elements.Nitrogen_3("N_3");

        private static final long serialVersionUID = 1L;
}
