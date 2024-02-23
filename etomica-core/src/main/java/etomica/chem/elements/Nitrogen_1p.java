package etomica.chem.elements;

    public class Nitrogen_1p extends ElementChemical {

        protected Nitrogen_1p(String symbol) {
            this(symbol, 14.01);
        }

        protected Nitrogen_1p(String symbol, double mass) {
            super(symbol, mass, 7);
        }

        public static final Nitrogen_1p INSTANCE = new Nitrogen_1p("N_1p");

        private static final long serialVersionUID = 1L;
}
