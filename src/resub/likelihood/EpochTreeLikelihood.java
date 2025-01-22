/*
* File TreeLikelihood.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/


package resub.likelihood;


import java.util.*;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.likelihood.BeerLikelihoodCore;
import beast.base.evolution.likelihood.BeerLikelihoodCore4;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;
import beast.base.inference.util.InputUtil;

@Description("Calculates the probability of sequence data on a beast.tree given a site and substitution model using " +
        "a variant of the 'peeling algorithm'. For details, see" +
        "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class EpochTreeLikelihood extends TreeLikelihood {
	public Input<List<SubstitutionModel>> modelsInput = new Input<>("model","substitution models that apply for certain time intervals", new ArrayList<SubstitutionModel>());
	public Input<Function> epochDatesInput = new Input<>("epochDates","list of threshold dates. " +
			"The list indicates the dates at which substitution models are switched.", Validate.REQUIRED);

	
	
    
	protected SubstitutionModel [] models;
	protected Function epochDates;
	protected double [] epochDatesArray;
	protected int nodeCount;
	// partialsMap maps node nr to partial that comes with it after taking epochs in account
	protected int [] partialsMap; 
	
	EpochBeagleTreeLikelihood epochbeagle;

    @Override
    public void initAndValidate() {
    	
        // sanity check: ensure the number of epoch dates is one less than the nr of models 
    	models = modelsInput.get().toArray(new SubstitutionModel[] {});
    	epochDates = epochDatesInput.get();
		if (models.length != epochDates.getDimension() + 1) {
			throw new IllegalArgumentException("The number of epoch dates should be one less than the number of substitution models");
		}
		epochDatesArray = new double[epochDates.getDimension()];
		
        // sanity check: alignment should have same #taxa as tree
		if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
			String leaves = "?";
			if (treeInput.get() instanceof Tree) {
				leaves = String.join(", ", ((Tree) treeInput.get()).getTaxaNames());
			}
			throw new IllegalArgumentException(String.format(
					"The number of leaves in the tree (%d) does not match the number of sequences (%d). "
							+ "The tree has leaves [%s], while the data refers to taxa [%s].",
					treeInput.get().getLeafNodeCount(), dataInput.get().getTaxonCount(),
					leaves, String.join(", ", dataInput.get().getTaxaNames())));
		}
		epochbeagle = null;
        
        
		epochbeagle = new EpochBeagleTreeLikelihood();
        try {
        	epochbeagle.initByName(
                    "data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(), 
                    "useTipLikelihoods", m_useTipLikelihoods.get(),"scaling", scaling.get().toString(),
                    "rootFrequencies", rootFrequenciesInput.get(),
                    "epochDates", epochDatesInput.get(),
                    "models", modelsInput.get());
	        if (epochbeagle.getBeagle() != null) {
	            //a Beagle instance was found, so we use it
	            return;
	        }
        } catch (Exception e) {
			// ignore
		}
        // No Beagle instance was found, so we use the good old java likelihood core
        epochbeagle = null;

        nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        // substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        int stateCount = dataInput.get().getMaxStateCount();
        int patterns = dataInput.get().getPatternCount();
        likelihoodCore = createLikelihoodCore(stateCount);

        String className = getClass().getSimpleName();

        Alignment alignment = dataInput.get();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if (proportionInvariant > 0) {
            calcConstantPatternIndices(patterns, stateCount);
        }

        initCore();

        patternLogLikelihoods = new double[patterns];
        m_fRootPartials = new double[patterns * stateCount];
        matrixSize = (stateCount + 1) * (stateCount + 1);
        probabilities = new double[(stateCount + 1) * (stateCount + 1)];
        Arrays.fill(probabilities, 1.0);

        if (dataInput.get().isAscertained) {
            useAscertainedSitePatterns = true;
        }
        
        partialsMap = new int[nodeCount];
    }

    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        int totalNodeCount = nodeCount * models.length + 1;
        likelihoodCore.initialize(
                totalNodeCount,
                dataInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1 + 1;
        final int intNodeCount = totalNodeCount - extNodeCount;

        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        } else {
            setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        }

        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i - 1);
        }

        int[] states = new int[dataInput.get().getPatternCount()];
       	Arrays.fill(states, dataInput.get().getDataType().getStateCount());
        likelihoodCore.setNodeStates(totalNodeCount - 1, states);
    }


    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;

    @Override
    public double calculateLogP() {
        if (epochbeagle != null) {
            logP = epochbeagle.calculateLogP();
            return logP;
        }
        final TreeInterface tree = treeInput.get();

        for (int i = 0;i < epochDatesArray.length; i++) {
        	epochDatesArray[i] = epochDates.getArrayValue(i);
        }
        for (int i = 0; i < nodeCount; i++) {
        	partialsMap[i] = i;
        }
        
        try {
        	if (traverse(tree.getRoot()) != Tree.IS_CLEAN)
        		calcLogP();
        }
        catch (ArithmeticException e) {
        	return Double.NEGATIVE_INFINITY;
        }
        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
            m_nScale = 0;
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            traverse(tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }

    /* Assumes there IS a branch rate model as opposed to traverse() */
    protected int traverse(final Node node) {

        int update = (node.isDirty() | hasDirt);

        int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        
        

		// First update the transition probability matrix(ices) for this branch
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
    		m_branchLengths[nodeIndex] = branchTime;

    		double startTime = node.getHeight();
    		double endTime = node.getParent().getHeight();
    		int startEpoch = Arrays.binarySearch(epochDatesArray, startTime);
    		if (startEpoch < 0) {
    			startEpoch = -startEpoch-1;
    		}
    		int endEpoch = Arrays.binarySearch(epochDatesArray, endTime);
    		if (endEpoch < 0) {
    			endEpoch = -endEpoch-1;
    		}

    		if (endEpoch != startEpoch) {
    			int childNum1 = nodeIndex;
        		for (int k = startEpoch; k < endEpoch; k++) {
                    likelihoodCore.setNodeMatrixForUpdate(childNum1);
                    for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                        final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                        models[k].getTransitionProbabilities(node, epochDatesArray[k], startTime, jointBranchRate, probabilities);
                        //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                        likelihoodCore.setNodeMatrix(childNum1, i, probabilities);
                    }    			
                    startTime = epochDatesArray[k];
                	likelihoodCore.setNodePartialsForUpdate(nodeIndex + (k+1) * nodeCount);
                	reportP(nodeIndex + (k+1) * nodeCount);
                    report(childNum1, (epochDatesArray.length+1)*nodeCount, nodeIndex + (k+1) * nodeCount);
                    likelihoodCore.calculatePartials(childNum1, (epochDatesArray.length+1)*nodeCount, nodeIndex + (k+1) * nodeCount);
                    childNum1 = nodeIndex + (k+1) * nodeCount;
        		}
        		nodeIndex = childNum1; 
    		}
    		partialsMap[node.getNr()] = nodeIndex;
    		
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                models[endEpoch].getTransitionProbabilities(node, endTime, startTime, jointBranchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
            	reportP(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculatePartials(partialsMap[childNum1], partialsMap[childNum2], nodeIndex);
                    report(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    likelihoodCore.integratePartials(nodeIndex, proportions, m_fRootPartials);

                    if (getConstantPattern() != null) { // && !SiteModel.g_bUseOriginal) {
                        proportionInvariant = m_siteModel.getProportionInvariant();
                        // some portion of sites is invariant, so adjust root partials for this
                        for (final int i : getConstantPattern()) {
                            m_fRootPartials[i] += proportionInvariant;
                        }
                    }

                    double[] rootFrequencies = models[models.length-1].getFrequencies();
                    if (rootFrequenciesInput.get() != null) {
                        rootFrequencies = rootFrequenciesInput.get().getFreqs();
                    }
                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
                }

            }
        }
        return update;
    } // traverse

    
    private void reportP(int i) {
//		System.err.println("flip " + i); 	
    }
    
    private void report(int childNum1, int childNum2, int nodeIndex) {
//		System.err.println(
//				(childNum1 < nodeCount ? partialsMap[childNum1]:childNum1) + " " + 
//				(childNum2 < nodeCount ? partialsMap[childNum2]:"") + " => " + 
//				nodeIndex);
	}

	@Override
    protected boolean requiresRecalculation() {
        if (epochbeagle != null) {
            return epochbeagle.requiresRecalculation();
        }
        
        
        if (InputUtil.isDirty(epochDatesInput)) {
        	hasDirt |= Tree.IS_DIRTY;
        	return super.requiresRecalculation() || true;
        }
        
        return super.requiresRecalculation();
    }

    @Override
    public void store() {
        if (epochbeagle != null) {
            epochbeagle.store();
            super.store();
            return;
        }
        super.store();
    }

    @Override
    public void restore() {
        if (epochbeagle != null) {
            epochbeagle.restore();
            super.restore();
            return;
        }
        super.restore();
    }

    
} // class TreeLikelihood
