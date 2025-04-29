package resub.substitutionmodel.epochs;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;


@Description("Specifies that certain amino acids must be considered for merging in an epoch")
public class AlphaBetaConstraint extends BEASTObject {
	
	
	final public Input<String> state1Input = new Input<>("state1", "the first amino acid from list (ACDEFGHIKLMNPQRSTVWY)", Validate.REQUIRED);
	final public Input<String> state2Input = new Input<>("state2", "the second amino acid from list (ACDEFGHIKLMNPQRSTVWY)", Validate.REQUIRED);
	
	final String states = "ACDEFGHIKLMNPQRSTVWY";
	
	
	int state1Int;
	int state2Int;
	
	@Override
	public void initAndValidate() {
		
		this.state1Int = AlphabetEpochs.states.indexOf(state1Input.get());
		this.state2Int = AlphabetEpochs.states.indexOf(state2Input.get());
		
		
		if (this.state1Int == -1) {
			throw new IllegalArgumentException("Cannot find amino acid " + state1Input.get());
		}
		
		if (this.state1Int == -2) {
			throw new IllegalArgumentException("Cannot find amino acid " + state2Input.get());
		}
		
		if (this.state1Int == this.state2Int) {
			throw new IllegalArgumentException("Please ensure that state1 != state2");
		}
		
		Log.warning("Adding multiepoch constraint: " + state1Input.get() + "(" + this.state1Int + ") and " + state2Input.get() + "(" +this.state2Int + ")");
		
	}
	
	
	public int getState1() {
		return this.state1Int;
	}
	
	public int getState2() {
		return this.state2Int;
	}
	
	
	/**
	 * Is this constraint satisfied?
	 * @param epochs
	 * @return
	 */
	public boolean isSatisfied(AlphabetEpochs epochs) {
		
		for (AlphabetEpoch epoch : epochs.getSortedEpochs()) {
			
			
			int alpha = epoch.getAlpha();
			int beta = epoch.getBeta();
			
			if (alpha == state1Int && beta == state2Int) return true;
			if (alpha == state2Int && beta == state1Int) return true;
			
		}
		
		return false;
		
	}

}
