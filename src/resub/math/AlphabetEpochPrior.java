package resub.math;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import resub.substitutionmodel.epochs.AlphabetEpoch;
import resub.substitutionmodel.epochs.AlphabetEpochs;


@Description("A prior on the epoch boundary parameters, mostly range validation checks")
public class AlphabetEpochPrior extends Distribution {
	
	final public Input<AlphabetEpochs> epochsInput = new Input<>("epochs", "vector of alphabet epochs", Input.Validate.REQUIRED);
	final public Input<RealParameter> pactiveInput = new Input<>("pactive", "proportion of epochs that follow this process", Input.Validate.REQUIRED);
	final public Input<RealParameter> prefineInput = new Input<>("prefine", "proportion of active epochs that are refinements", Input.Validate.REQUIRED);

	final public Input<Boolean> checkStatesInput = new Input<>("checkStates", "turn off to permit any combination of states (debugging)", false);

	
	
	AlphabetEpochs epochs;
	
	@Override
    public void initAndValidate() {
		
		this.epochs = epochsInput.get();
        calculateLogP();
    }

    @Override
    public double calculateLogP() {
    	
    	
    	logP = 0;
    	

    	
    	
    	if (checkStatesInput.get()) {
	    
    		// Get sorted epochs
        	List<AlphabetEpoch> epochList = epochs.getSortedEpochs();
    		
	    	// Something has gone wrong
	    	if (epochList == null) {
	    		//Log.warning("A");
	    		logP = Double.NEGATIVE_INFINITY;
				return logP;
	    	}
	    	
	    	
	    	// Check the indicators are in range
	    	for(int i = 0; i < epochList.size(); i ++) {
	    		AlphabetEpoch epoch = epochList.get(i);
	    		if (epoch.getIndicator() < 0 || epoch.getIndicator() > 2){
	    			//Log.warning("B");
	    			logP = Double.NEGATIVE_INFINITY;
	    			return logP;
	    		}
	    	}
	    	
	    	
	    	// Ensure that all amino acid constraints are satisfied
	    	if (!epochs.constraintsAreSatisfied()) {
	    		logP = Double.NEGATIVE_INFINITY;
    			return logP;
	    	}
	    	
	    	for(int i = 0; i < epochList.size(); i ++) {
	    		
	    		
	    		AlphabetEpoch epoch = epochList.get(i);
	    		
	    		// Cannot merge 2 states if they are the same
	    		if (epoch.getAlpha() == epoch.getBeta()) {
	    			//Log.warning("C");
	    			logP = Double.NEGATIVE_INFINITY;
	    			return logP;
	    		}
	    		
	    			
				// Check that alpha is not pointing to an invalid state
				String stateNameInEpoch = epoch.getStateNames().get(epoch.getAlpha());
				if (stateNameInEpoch.equals(AlphabetEpoch.BLANK_STATE)) {
					logP = Double.NEGATIVE_INFINITY;
					//Log.warning("Reject: alpha is pointing to a deleted state");
	    			return logP;
				}
				
				
				// The probability of seeing alpha,beta is 1 over the total number of valid states, conditional on the alphabet size at this epoch
				int s = epoch.getNStatesAtEpoch();
				logP -= Math.log(s * (s-1));
	    		
	    	}
	    	
	    	
    	}
    	
    	
    	// Refinement expansion probabilities
    	double pactive = pactiveInput.get().getValue();
    	//double prefine = pactive * (prefineInput.get().getValue() / 2); // Halved to account for the fact that (a,b) is the same as (b,a)
    	double prefine = pactive * prefineInput.get().getValue(); 
    	double pexpand = (pactive - prefine);
    	
    	if (pactive < 0 || pactive > 1 || prefine < 0 || prefine > 1 || pexpand < 0 || pexpand > 1) {
    		//Log.warning("E");
    		logP = Double.NEGATIVE_INFINITY;
			return logP;
    	}
    	
    	int n = epochs.getNOldEpochs();
    	for(int i = 0; i < n; i ++) {
    	
    		AlphabetEpoch epoch = epochs.getEpochByIndex(i);
    		if (epoch.getIndicator() == 0) {
    			if (pactive == 1) logP = Double.NEGATIVE_INFINITY;
    			logP += Math.log(1-pactive);
    		}
    		
    		else if (epoch.getIndicator() == 1) {
    			if (prefine == 0) logP = Double.NEGATIVE_INFINITY;
    			logP += Math.log(prefine); 
    		}
    		
    		else if (epoch.getIndicator() == 2 ) {
    			if (pexpand == 0) logP = Double.NEGATIVE_INFINITY;
    			logP += Math.log(pexpand);
    		}
    		
    		
    	}
    	
    	
    	//Log.warning("F" + logP);
    	return logP;
    	
    }
	
	@Override
	public List<String> getArguments() {
		List<String> args = new ArrayList<>();
		args.add(this.epochs.getAlphas().getID());
		args.add(this.epochs.getBetas().getID());
		args.add(this.epochs.getIndicators().getID());
		return args;
	}

	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		conds.add(this.epochs.getTransitionAges().getID());
		conds.add(pactiveInput.get().getID());
		conds.add(prefineInput.get().getID());
		return conds;
	}

	@Override
	public void sample(State state, Random random) {
		
		if (sampledFlag) return;
        sampledFlag = true;
        
        System.out.println("Sampling");
        
        // Sample conditions
        sampleConditions(state, random);
        
        IntegerParameter indicators = epochs.getIndicators();
        double pactive = pactiveInput.get().getArrayValue();
        //double prefineGivenActive = prefineInput.get().getValue() / 2; // Halved to account for the fact that (a,b) is the same as (b,a)
        double prefineGivenActive = prefineInput.get().getValue(); 
        
        
        // Sample the indicators
    	for(int i = 0; i < epochs.getNOldEpochs(); i ++) {
    		
    		// Active
    		if (random.nextDouble() < pactive) {
    			
    			// Refinement
    			if (random.nextDouble() < prefineGivenActive) {
    				indicators.setValue(i, 1);
    			}
    			
    			// Expansion
    			else {
    				indicators.setValue(i, 2);
    			}
    		
    		}
    		
    		// Inactive
    		else {
    			indicators.setValue(i, 0);
    		}
    		
		}
    	
//    	System.out.print("I = ");
//    	for(int i = 0; i < epochs.getNOldEpochs(); i ++) {
//    		System.out.print(epochs.getEpochByIndex(i).getIndicator() + " ");
//    	}
//    	System.out.println();
//    	

    	
    	// Sample alpha and beta, going from bottom to top
    	int nstates = epochs.getNStatesTotal();
    	List<String> stateList = epochs.getStateList();
    	
    	for(int i = 0; i < epochs.getNOldEpochs(); i ++) {
    		
    		epochs.requestUpdate();
    		AlphabetEpoch epoch = epochs.getEpochByRank(i);
    		
			// Sample alpha and beta such that alpha!=beta, and that neither index corresponds to a blank state in this epoch
			int alpha = sampleState(nstates, stateList, random, -1);
			int beta = sampleState(nstates, stateList, random, alpha);
			
			/*
			System.out.println("rank " + i + "; alpha=" + alpha + "; beta=" + beta);
			System.out.println("indicator=" + epoch.getIndicator());
			System.out.println("age=" + epoch.getAge());
			System.out.println("alpha=" + stateList.get(alpha));
			System.out.println("beta=" + stateList.get(beta));
			*/

			
			epoch.setAlpha(alpha);
			epoch.setBeta(beta);
			
			// Update state list
			String alphaState = stateList.get(alpha);
			String betaState = stateList.get(beta);
			if (epoch.isRefinement()) {
				String mergedStateName = alphaState + "/" + betaState;
				stateList.set(stateList.indexOf(alphaState), mergedStateName);
				stateList.set(stateList.indexOf(betaState), ".");
			}
			else if (epoch.isExpansion()) {
				stateList.set(stateList.indexOf(betaState), ".");
			}
    			
    		
    		
//			System.out.print("States: ");
//			for (String str : stateList) {
//				System.out.print(str + ",");
//			}
//			System.out.println();
    		
    	}
    	
//    	
//    	System.out.print("alpha = ");
//    	for(int i = 0; i < epochs.getNOldEpochs(); i ++) {
//    		System.out.print(epochs.getEpochByIndex(i).getAlpha() + " ");
//    	}
//    	System.out.println();
//    	
//    	System.out.print("beta = ");
//    	for(int i = 0; i < epochs.getNOldEpochs(); i ++) {
//    		System.out.print(epochs.getEpochByIndex(i).getBeta() + " ");
//    	}
//    	System.out.println();
        
		
	}
	
	
	// Keep sampling a state until it is not a blank state (.)
	private int sampleState(int nstates, List<String> possibleStates, Random random, int taboo) {
		int val = random.nextInt(nstates);
		while (val == taboo || possibleStates.get(val).equals(AlphabetEpoch.BLANK_STATE)){
			val = random.nextInt(nstates);
		}
		return val;
	}

}
