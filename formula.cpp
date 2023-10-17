// This file contains code from NVSim, (c) 2012-2013,  Pennsylvania State University 
//and Hewlett-Packard Company. See LICENSE_NVSim file in the top-level directory.
//No part of DESTINY Project, including this file, may be copied,
//modified, propagated, or distributed except according to the terms
//contained in the LICENSE file.

#include "formula.h"
#include "constant.h"
#include <stdlib.h>

bool isPow2(int n) {
	if (n < 1)
		return false;
	return !(n & (n - 1));
}

double CalculateGateCap(double width, Technology tech) {
	return (tech.capIdealGate + tech.capOverlap + 3 * tech.capFringe) * width
			+ tech.phyGateLength * tech.capPolywire;
}

double CalculateFBRAMGateCap(double width, double thicknessFactor, Technology tech) {
	return (tech.capIdealGate / thicknessFactor + tech.capOverlap + 3 * tech.capFringe) * width
			+ tech.phyGateLength * tech.capPolywire;
}

double CalculateFBRAMDrainCap(double width, Technology tech) {
	return (3 * tech.capSidewall + tech.capDrainToChannel) * width;
}

double CalculateGateArea(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double heightTransistorRegion, Technology tech,
		double *height, double *width) {
	double	ratio = widthPMOS / (widthPMOS + widthNMOS);

	double maxWidthPMOS, maxWidthNMOS;
	double unitWidthRegionP, unitWidthRegionN;
	double widthRegionP, widthRegionN;
	double heightRegionP, heightRegionN;

	if (ratio == 0) {	/* no PMOS */
		maxWidthPMOS = 0;
		maxWidthNMOS = heightTransistorRegion;
	} else if (ratio == 1) {	/* no NMOS */
		maxWidthPMOS = heightTransistorRegion;
		maxWidthNMOS = 0;
	} else {
		maxWidthPMOS = ratio * (heightTransistorRegion - MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize);
		maxWidthNMOS = maxWidthPMOS / ratio * (1 - ratio);
	}

	if (widthPMOS > 0) {
		if (widthPMOS < maxWidthPMOS) { /* No folding */
			unitWidthRegionP = tech.featureSize;
			heightRegionP = widthPMOS;
		} else {	/* Folding */
			int numFoldedPMOS = (int)(ceil(widthPMOS / (maxWidthPMOS - 3 * tech.featureSize)));	/* 3F for folding overhead */
			unitWidthRegionP = numFoldedPMOS * tech.featureSize + (numFoldedPMOS-1) * tech.featureSize * MIN_GAP_BET_POLY;
			heightRegionP = maxWidthPMOS;
		}
	} else {
		unitWidthRegionP = 0;
		heightRegionP = 0;
	}

	if (widthNMOS > 0) {
		if (widthNMOS < maxWidthNMOS) { /* No folding */
			unitWidthRegionN = tech.featureSize;
			heightRegionN = widthNMOS;
		} else {	/* Folding */
			int numFoldedNMOS = (int)(ceil(widthNMOS / (maxWidthNMOS - 3 * tech.featureSize)));	/* 3F for folding overhead */
			unitWidthRegionN = numFoldedNMOS * tech.featureSize + (numFoldedNMOS-1) * tech.featureSize * MIN_GAP_BET_POLY;
			heightRegionN = maxWidthNMOS;
		}
	} else {
		unitWidthRegionN = 0;
		heightRegionN = 0;
	}

	switch (gateType) {
	case INV:
		widthRegionP = 2 * tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2) + unitWidthRegionP;
		widthRegionN = 2 * tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2) + unitWidthRegionN;
		break;
	case NOR:
		widthRegionP = 2 * tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2)
						+ unitWidthRegionP * numInput + (numInput - 1) * tech.featureSize * MIN_GAP_BET_POLY;
		widthRegionN = 2 * tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2)
						+ unitWidthRegionN * numInput
						+ (numInput - 1) * tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2);
		break;
	case NAND:
		widthRegionN = 2 * tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2)
						+ unitWidthRegionN * numInput + (numInput - 1) * tech.featureSize * MIN_GAP_BET_POLY;
		widthRegionP = 2 * tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2)
						+ unitWidthRegionP * numInput
						+ (numInput - 1) * tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2);
		break;
	default:
		widthRegionN = widthRegionP = 0;
	}

	*width = MAX(widthRegionN, widthRegionP);
	if (widthPMOS > 0 && widthNMOS > 0) {	/* it is a gate */
		*height = heightRegionN + heightRegionP + tech.featureSize * MIN_GAP_BET_P_AND_N_DIFFS
					+ 2 * tech.featureSize * MIN_WIDTH_POWER_RAIL;
	} else {	/* it is a transistor */
		*height = heightRegionN + heightRegionP;	/* one of them is zero, and no power rail is added */
	}

	return (*width)*(*height);
}

void CalculateGateCapacitance(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double heightTransistorRegion, Technology tech,
		double *capInput, double *capOutput) {

	if(tech.featureSize >= 22){
	/* TO-DO: most parts of this function is the same of CalculateGateArea,
	 * perhaps they will be combined in future
	 */
		double	ratio = widthPMOS / (widthPMOS + widthNMOS);

		double maxWidthPMOS = 0, maxWidthNMOS = 0;
		double unitWidthDrainP = 0, unitWidthDrainN = 0;
		double widthDrainP = 0, widthDrainN = 0;
		double heightDrainP = 0, heightDrainN = 0;
		int numFoldedPMOS = 1, numFoldedNMOS = 1;

		if (ratio == 0) {	/* no PMOS */
			maxWidthPMOS = 0;
			maxWidthNMOS = heightTransistorRegion;
		} else if (ratio == 1) {	/* no NMOS */
			maxWidthPMOS = heightTransistorRegion;
			maxWidthNMOS = 0;
		} else {
			maxWidthPMOS = ratio * (heightTransistorRegion - MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize);
			maxWidthNMOS = maxWidthPMOS / ratio * (1 - ratio);
		}

		if (widthPMOS > 0) {
			if (widthPMOS < maxWidthPMOS) { /* No folding */
				unitWidthDrainP = 0;
				heightDrainP = widthPMOS;
			} else {	/* Folding */
				if (maxWidthPMOS < 3 * tech.featureSize) {
					cout << "Error: Unable to do PMOS folding because PMOS size limitation is less than 3F!" <<endl;
					exit(-1);
				}
				numFoldedPMOS = (int)(ceil(widthPMOS / (maxWidthPMOS - 3 * tech.featureSize)));	/* 3F for folding overhead */
				unitWidthDrainP = (numFoldedPMOS-1) * tech.featureSize * MIN_GAP_BET_POLY;
				heightDrainP = maxWidthPMOS;
			}
		} else {
			unitWidthDrainP = 0;
			heightDrainP = 0;
		}

		if (widthNMOS > 0) {
			if (widthNMOS < maxWidthNMOS) { /* No folding */
				unitWidthDrainN = 0;
				heightDrainN = widthNMOS;
			} else {	/* Folding */
				if (maxWidthNMOS < 3 * tech.featureSize) {
					cout << "Error: Unable to do NMOS folding because NMOS size limitation is less than 3F!" <<endl;
					exit(-1);
				}
				numFoldedNMOS = (int)(ceil(widthNMOS / (maxWidthNMOS - 3 * tech.featureSize)));	/* 3F for folding overhead */
				unitWidthDrainN = (numFoldedNMOS-1) * tech.featureSize * MIN_GAP_BET_POLY;
				heightDrainN = maxWidthNMOS;
			}
		} else {
			unitWidthDrainN = 0;
			heightDrainN = 0;
		}

		switch (gateType) {
		case INV:
			if (widthPMOS > 0)
				widthDrainP = tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2) + unitWidthDrainP;
			if (widthNMOS > 0)
				widthDrainN = tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2) + unitWidthDrainN;
			break;
		case NOR:
			/* PMOS is in series, worst case capacitance is below */
			if (widthPMOS > 0)
				widthDrainP = tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2)
							+ unitWidthDrainP * numInput + (numInput - 1) * tech.featureSize * MIN_GAP_BET_POLY;
			/* NMOS is parallel, capacitance is multiplied as below */
			if (widthNMOS > 0)
				widthDrainN = (tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2)
							+ unitWidthDrainN) * numInput;
			break;
		case NAND:
			/* NMOS is in series, worst case capacitance is below */
			if (widthNMOS > 0)
				widthDrainN = tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2)
							+ unitWidthDrainN * numInput + (numInput - 1) * tech.featureSize * MIN_GAP_BET_POLY;
			/* PMOS is parallel, capacitance is multiplied as below */
			if (widthPMOS > 0)
				widthDrainP = (tech.featureSize * (CONTACT_SIZE + MIN_GAP_BET_CONTACT_POLY * 2)
							+ unitWidthDrainP) * numInput;
			break;
		default:
			widthDrainN = widthDrainP = 0;
		}

		/* Junction capacitance */
		double capDrainBottomN = widthDrainN * heightDrainN * tech.capJunction;
		double capDrainBottomP = widthDrainP * heightDrainP * tech.capJunction;

		/* Sidewall capacitance */
		double capDrainSidewallN, capDrainSidewallP;
		if (numFoldedNMOS % 2 == 0)
			capDrainSidewallN = 2 * widthDrainN * tech.capSidewall;
		else
			capDrainSidewallN = (2 * widthDrainN + heightDrainN) * tech.capSidewall;
		if (numFoldedPMOS % 2 == 0)
			capDrainSidewallP = 2 * widthDrainP * tech.capSidewall;
		else
			capDrainSidewallP = (2* widthDrainP + heightDrainP) * tech.capSidewall;

		/* Drain to channel capacitance */
		double capDrainToChannelN = numFoldedNMOS * heightDrainN * tech.capDrainToChannel;
		double capDrainToChannelP = numFoldedPMOS * heightDrainP * tech.capDrainToChannel;

		if (capOutput)
			*(capOutput) = capDrainBottomN + capDrainBottomP + capDrainSidewallN + capDrainSidewallP + capDrainToChannelN + capDrainToChannelP;
		if (capInput)
			*(capInput) = CalculateGateCap(widthNMOS, tech) + CalculateGateCap(widthPMOS, tech);
	} else {
		int speciallayout = 0;
		if (heightTransistorRegion != tech.featureSize *MAX_TRANSISTOR_HEIGHT) speciallayout = 1;


		if (capInput){
			*(capInput) = CalculateGateCap(widthNMOS, tech) + CalculateGateCap(widthPMOS, tech);
		}
		if (capOutput){
			if (tech.featureSize <= 14) {  // 1.4 update: finfet or GAA

				if (tech.featureSize > 2) { // 1.4 update: finfet
					widthNMOS *= 1/(2 * tech.featureSize); 
					widthPMOS *= 1/(2 * tech.featureSize); 
				}

				else if (tech.featureSize <= 2) { // 1.4 update: GAA 
					widthNMOS *= 1/(2 * tech.featureSize);
					widthPMOS *= 1/(2 * tech.featureSize);
				}

				// 1.4 update: handle updated standard cell trend below 14 nm node

				heightTransistorRegion *= ((double) MAX_TRANSISTOR_HEIGHT_FINFET/MAX_TRANSISTOR_HEIGHT);

				if (tech.featureSize == 14)
				heightTransistorRegion *= ( (double)MAX_TRANSISTOR_HEIGHT_14nm /MAX_TRANSISTOR_HEIGHT_FINFET);
				else if (tech.featureSize == 10)
				heightTransistorRegion *= ( (double)MAX_TRANSISTOR_HEIGHT_10nm /MAX_TRANSISTOR_HEIGHT_FINFET);
				else if (tech.featureSize == 7)
				heightTransistorRegion *= ( (double)MAX_TRANSISTOR_HEIGHT_7nm /MAX_TRANSISTOR_HEIGHT_FINFET);
				else if (tech.featureSize == 5)
				heightTransistorRegion *= ( (double)MAX_TRANSISTOR_HEIGHT_5nm /MAX_TRANSISTOR_HEIGHT_FINFET);
				else if (tech.featureSize == 3)
				heightTransistorRegion *= ( (double)MAX_TRANSISTOR_HEIGHT_3nm /MAX_TRANSISTOR_HEIGHT_FINFET);
				else if (tech.featureSize == 2)
				heightTransistorRegion *= ( (double)MAX_TRANSISTOR_HEIGHT_2nm /MAX_TRANSISTOR_HEIGHT_FINFET);
				else if (tech.featureSize == 1)
				heightTransistorRegion *= ( (double)MAX_TRANSISTOR_HEIGHT_1nm /MAX_TRANSISTOR_HEIGHT_FINFET);
				else
				heightTransistorRegion *= 1;
			}

			double	ratio = widthPMOS / (widthPMOS + widthNMOS);
			double maxWidthPMOS = 0, maxWidthNMOS = 0;
			int maxNumPFin = 0, maxNumNFin = 0;	/* Max numbers of fin for the specified cell height */
			double unitWidthDrainP = 0, unitWidthDrainN = 0;
			double unitWidthSourceP = 0, unitWidthSourceN = 0;
			double widthDrainP = 0, widthDrainN = 0;
			double heightDrainP = 0, heightDrainN = 0;
			int numFoldedPMOS = 1, numFoldedNMOS = 1;
			double widthDrainSidewallP = 0, widthDrainSidewallN = 0;
			
			
			// 1.4 update: add GAA-related parameters
			int NumPSheet;
			int NumNSheet;
			int maxNumSheet;
			int maxNumPSheet; 
			int maxNumNSheet; 

			if (tech.featureSize >= 22) { // Bulk
				if (ratio == 0) {	/* no PMOS */
					maxWidthPMOS = 0;
					maxWidthNMOS = heightTransistorRegion - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize;
				} else if (ratio == 1) {	/* no NMOS */
					maxWidthPMOS = heightTransistorRegion - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize;
					maxWidthNMOS = 0;
				} else {
					maxWidthPMOS = ratio * (heightTransistorRegion - MIN_GAP_BET_P_AND_N_DIFFS * tech.featureSize - (MIN_POLY_EXT_DIFF + MIN_GAP_BET_FIELD_POLY/2) * 2 * tech.featureSize);
					maxWidthNMOS = maxWidthPMOS / ratio * (1 - ratio);
				}

				if (widthPMOS > 0) {
					if (widthPMOS <= maxWidthPMOS) { /* No folding */
						unitWidthDrainP = tech.featureSize * MIN_GAP_BET_GATE_POLY;
						unitWidthSourceP = unitWidthDrainP;
						heightDrainP = widthPMOS;
					} else {	/* Folding */
						numFoldedPMOS = (int)(ceil(widthPMOS / maxWidthPMOS));
						unitWidthDrainP = (int)ceil((double)(numFoldedPMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;	// Num of drain fingers >= num of source fingers
						unitWidthSourceP = (int)floor((double)(numFoldedPMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;
						heightDrainP = maxWidthPMOS;
					}
				} else {
					unitWidthDrainP = 0;
					unitWidthSourceP = 0;
					heightDrainP = 0;
				}
				if (widthNMOS > 0) {
					if (widthNMOS <= maxWidthNMOS) { /* No folding */
						unitWidthDrainN = tech.featureSize * MIN_GAP_BET_GATE_POLY;
						unitWidthSourceN = unitWidthDrainN;
						heightDrainN = widthNMOS;
					} else {	/* Folding */
						numFoldedNMOS = (int)(ceil(widthNMOS / maxWidthNMOS));
						unitWidthDrainN = (int)ceil((double)(numFoldedNMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;
						unitWidthSourceN = (int)floor((double)(numFoldedNMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY;
						heightDrainN = maxWidthNMOS;
					}
				} else {
					unitWidthDrainN = 0;
					unitWidthSourceN = 0;
					heightDrainN = 0;
				}

			} else { // 1.4 update: FinFET or GAA

				// 1.4 update: add variables related to cell width
				double modified_POLY_WIDTH_FINFET= POLY_WIDTH_FINFET;
				double CPP_advanced = POLY_WIDTH_FINFET + MIN_GAP_BET_GATE_POLY_FINFET;
				double modified_MIN_GAP_BET_GATE_POLY_FINFET;

				// 1.4 update: handle different fin number/CPP/Cell width trend for 14 nm and below
				if (tech.featureSize == 14) { // adding more cases 
					maxNumPFin = maxNumNFin = tech.max_fin_num; // changed from 3 to 4
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_14nm;
					CPP_advanced = CPP_14nm;
				} else if (tech.featureSize == 10) {
					maxNumPFin = maxNumNFin = tech.max_fin_num;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_10nm;
					CPP_advanced = CPP_10nm;
				} else if (tech.featureSize == 7) {
					maxNumPFin = maxNumNFin = tech.max_fin_num;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_7nm;
					CPP_advanced = CPP_7nm;
				} else if (tech.featureSize == 5) {
					maxNumPFin = maxNumNFin = tech.max_fin_num;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_5nm;
					CPP_advanced = CPP_5nm;
				} else if (tech.featureSize == 3) {
					maxNumPFin = maxNumNFin = tech.max_fin_num;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_3nm;
					CPP_advanced = CPP_3nm;
				} else if (tech.featureSize == 2) {
					maxNumPSheet= maxNumNSheet = tech.max_fin_per_GAA;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_2nm;
					CPP_advanced = CPP_2nm;
				} else if (tech.featureSize == 1) {
					maxNumPSheet= maxNumNSheet = tech.max_fin_per_GAA;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_1nm;
					CPP_advanced = CPP_1nm;
				} 


			// 1.4 update: setting the maximum number of fins
			if (!speciallayout) {

				if (tech.featureSize == 14) { // adding more cases 
					maxNumPFin = maxNumNFin = tech.max_fin_num; // changed from 3 to 4
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_14nm;
					CPP_advanced = CPP_14nm;
				} else if (tech.featureSize == 10) {
					maxNumPFin = maxNumNFin = tech.max_fin_num;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_10nm;
					CPP_advanced = CPP_10nm;
				} else if (tech.featureSize == 7) {
					maxNumPFin = maxNumNFin = tech.max_fin_num;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_7nm;
					CPP_advanced = CPP_7nm;
				} else if (tech.featureSize == 5) {
					maxNumPFin = maxNumNFin = tech.max_fin_num;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_5nm;
					CPP_advanced = CPP_5nm;
				} else if (tech.featureSize == 3) {
					maxNumPFin = maxNumNFin = tech.max_fin_num;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_3nm;
					CPP_advanced = CPP_3nm;
				} else if (tech.featureSize == 2) {
					maxNumPSheet= maxNumNSheet = tech.max_fin_per_GAA;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_2nm;
					CPP_advanced = CPP_2nm;
				} else if (tech.featureSize == 1) {
					maxNumPSheet= maxNumNSheet = tech.max_fin_per_GAA;
					modified_POLY_WIDTH_FINFET= POLY_WIDTH_1nm;
					CPP_advanced = CPP_1nm;
				} 
			}

			else {

			if (tech.featureSize == 14) { 
				maxNumPFin = maxNumNFin =  (floor)( ( (heightTransistorRegion - (double) MIN_GAP_BET_P_AND_N_DIFFS_14nm * tech.featureSize -  (double) OUTER_HEIGHT_REGION_14nm * tech.featureSize) + tech.PitchFin ) / 2.0 / (tech.widthFin + tech.PitchFin) );
				modified_POLY_WIDTH_FINFET= POLY_WIDTH_14nm;
				CPP_advanced = CPP_14nm;
		
			} else if (tech.featureSize == 10) {
				maxNumPFin = maxNumNFin =  (floor)( ( (heightTransistorRegion -  (double) MIN_GAP_BET_P_AND_N_DIFFS_10nm * tech.featureSize -  (double) OUTER_HEIGHT_REGION_10nm * tech.featureSize) + tech.PitchFin ) / 2.0 / (tech.widthFin + tech.PitchFin) );
				modified_POLY_WIDTH_FINFET= POLY_WIDTH_10nm;
				CPP_advanced = CPP_10nm;        
			
			} else if (tech.featureSize == 7) {
				maxNumPFin = maxNumNFin = (floor)( ( (heightTransistorRegion -  (double) MIN_GAP_BET_P_AND_N_DIFFS_7nm * tech.featureSize -  (double) OUTER_HEIGHT_REGION_7nm * tech.featureSize) + tech.PitchFin ) / 2.0 / (tech.widthFin + tech.PitchFin) );
				modified_POLY_WIDTH_FINFET= POLY_WIDTH_7nm;
				CPP_advanced = CPP_7nm;         
			
			} else if (tech.featureSize == 5) {
				maxNumPFin = maxNumNFin = (floor)( ( (heightTransistorRegion -  (double) MIN_GAP_BET_P_AND_N_DIFFS_5nm * tech.featureSize -  (double) OUTER_HEIGHT_REGION_5nm * tech.featureSize) + tech.PitchFin ) / 2.0 / (tech.widthFin + tech.PitchFin) );
				modified_POLY_WIDTH_FINFET= POLY_WIDTH_5nm;
				CPP_advanced = CPP_5nm;         
			
			} else if (tech.featureSize == 3) {
				maxNumPFin = maxNumNFin = (floor)( ( (heightTransistorRegion -  (double) MIN_GAP_BET_P_AND_N_DIFFS_3nm * tech.featureSize -  (double) OUTER_HEIGHT_REGION_3nm * tech.featureSize) + tech.PitchFin ) / 2 / (tech.widthFin + tech.PitchFin) );
				modified_POLY_WIDTH_FINFET= POLY_WIDTH_3nm;
				CPP_advanced = CPP_3nm;           
			
			} else if (tech.featureSize == 2) {
				maxNumPSheet= maxNumNSheet = (floor)( ( (heightTransistorRegion -  (double) MIN_GAP_BET_P_AND_N_DIFFS_2nm * tech.featureSize -  (double) OUTER_HEIGHT_REGION_2nm* tech.featureSize) + tech.PitchFin ) / 2 / (tech.widthFin + tech.PitchFin) );
				modified_POLY_WIDTH_FINFET= POLY_WIDTH_2nm;
				CPP_advanced = CPP_2nm;           
			
			} else if (tech.featureSize == 1) {
				maxNumPSheet= maxNumNSheet = (floor)( ( (heightTransistorRegion -  (double) MIN_GAP_BET_P_AND_N_DIFFS_1nm * tech.featureSize -  (double) OUTER_HEIGHT_REGION_1nm * tech.featureSize) + tech.PitchFin ) / 2 / (tech.widthFin + tech.PitchFin) );
				modified_POLY_WIDTH_FINFET= POLY_WIDTH_1nm;
				CPP_advanced = CPP_1nm;           
			
			} 
				

			}

			// 1.4 update: temp_P, temp_N, temp_P_NS, temp_N_NS
			// 1.4 update: temp_N_ratio, temp_P_ratio, temp_N_NS_ratio, temp_P_NS_ratio

			double temp_P=2*maxNumPFin;
			double temp_N=2*maxNumNFin;
			double temp_P_NS=2*maxNumPSheet;
			double temp_N_NS=2*maxNumNSheet;

			int temp_N_ratio;
			int temp_P_ratio;
			int temp_N_NS_ratio;
			int temp_P_NS_ratio;

			// 1.4 update: setting the maximum number of fin for folding

			if (ratio == 0) {   /* no PFinFET */

				maxNumPFin = 0;
				maxNumNFin = temp_N;
				maxNumPSheet = 0;
				maxNumNSheet= temp_N_NS;

			} else if (ratio == 1) {    /* no NFinFET */

				maxNumPFin = temp_P;
				maxNumNFin = 0;
				maxNumPSheet = temp_P_NS;
				maxNumNSheet= 0;

			} else {

				if (ratio>0.5) {
				temp_N_ratio=int (temp_N*(1-ratio));
				temp_P_ratio= temp_N-temp_N_ratio;
				temp_N_NS_ratio=tech.max_fin_per_GAA; // maximum fin for GAA is fixed to 1 
				temp_P_NS_ratio= tech.max_fin_per_GAA; // maximum fin for GAA is fixed to 1
				}

				else {
				temp_P_ratio=int (temp_P*(ratio));
				temp_N_ratio= temp_P - temp_P_ratio;
				temp_P_NS_ratio=tech.max_fin_per_GAA; // maximum fin for GAA is fixed to 1
				temp_N_NS_ratio=tech.max_fin_per_GAA; // maximum fin for GAA is fixed to 1
				}

				if (temp_P_ratio==0) {temp_P_ratio +=1; temp_N_ratio = 2*maxNumPFin-temp_P_ratio;} // handle max fin number = 0 casees, since they don't make sense 
				if (temp_N_ratio==0) {temp_N_ratio +=1; temp_P_ratio = 2*maxNumNFin-temp_N_ratio;} // handle max fin number = 0 casees, since they don't make sense
				if (temp_P_NS_ratio==0) {temp_P_NS_ratio +=1; temp_N_NS_ratio = 2*maxNumPSheet-temp_P_NS_ratio;} // handle max fin number = 0 casees, since they don't make sense
				if (temp_N_NS_ratio==0) {temp_N_NS_ratio +=1; temp_P_NS_ratio = 2*maxNumNSheet-temp_N_NS_ratio;} // handle max fin number = 0 casees, since they don't make sense

				maxNumPFin = temp_P_ratio;
				maxNumNFin = temp_N_ratio;
				maxNumPSheet = temp_P_NS_ratio; 
				maxNumNSheet = temp_N_NS_ratio; 
			}

				// 1.4 update: updated the gap distance between finfets
				modified_MIN_GAP_BET_GATE_POLY_FINFET= CPP_advanced-modified_POLY_WIDTH_FINFET; 

				// 1.4 update: capacitance calculation for 14 nm and beyond
				if (tech.featureSize < 22  && tech.featureSize >= 3) { // 1.4 update: FinFET case

					int NumPFin = (int)(ceil(widthPMOS));
					if (NumPFin > 0) {
						if (NumPFin <= maxNumPFin) { /* No folding */
							unitWidthDrainP = tech.featureSize * modified_MIN_GAP_BET_GATE_POLY_FINFET;
							unitWidthSourceP = unitWidthDrainP;
							heightDrainP = NumPFin * tech.widthFin;
						} else {    /* Folding */
							numFoldedPMOS = (int)(ceil(NumPFin / maxNumPFin));
							unitWidthDrainP = (int)ceil((double)(numFoldedPMOS+1)/2) * tech.featureSize * modified_MIN_GAP_BET_GATE_POLY_FINFET;
							unitWidthSourceP = (int)floor((double)(numFoldedPMOS+1)/2) * tech.featureSize * modified_MIN_GAP_BET_GATE_POLY_FINFET;
							heightDrainP = maxNumPFin * tech.widthFin; // modified from heightDrainP = (maxNumPFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
						}
					} else {
						unitWidthDrainP = 0;
						unitWidthSourceP = 0;
						heightDrainP = 0;
					}

					int NumNFin = (int)(ceil(widthNMOS));

					if (NumNFin > 0) {
						if (NumNFin <= maxNumNFin) { /* No folding */
							unitWidthDrainN = tech.featureSize * modified_MIN_GAP_BET_GATE_POLY_FINFET;
							unitWidthSourceN = unitWidthDrainN;
							heightDrainN = NumNFin * tech.widthFin;
						} else {    /* Folding */
							numFoldedNMOS = (int)(ceil(NumNFin / maxNumNFin));
							unitWidthDrainN = (int)ceil((double)(numFoldedNMOS+1)/2) * tech.featureSize * modified_MIN_GAP_BET_GATE_POLY_FINFET;
							unitWidthSourceN = (int)floor((double)(numFoldedNMOS+1)/2) * tech.featureSize * modified_MIN_GAP_BET_GATE_POLY_FINFET;
							heightDrainN = maxNumNFin * tech.widthFin; // modified from heightDrainN = (maxNumNFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
						}
					} else {
						unitWidthDrainN = 0;
						unitWidthSourceN = 0;
						heightDrainN = 0;
					}
				}

				else {
					int NumPSheet = (int)(ceil(widthPMOS));

					if (NumPSheet > 0) {
						if (NumPSheet <= maxNumPSheet) { /* No folding */
							unitWidthDrainP = tech.featureSize * modified_MIN_GAP_BET_GATE_POLY_FINFET;
							unitWidthSourceP = unitWidthDrainP;
							heightDrainP = NumPSheet * tech.widthFin;
						} else {    /* Folding */
							numFoldedPMOS = (int)(ceil(NumPSheet/ maxNumPSheet));
							unitWidthDrainP = (int)ceil((double)(numFoldedPMOS+1)/2) * tech.featureSize * modified_MIN_GAP_BET_GATE_POLY_FINFET;
							unitWidthSourceP = (int)floor((double)(numFoldedPMOS+1)/2) * tech.featureSize * modified_MIN_GAP_BET_GATE_POLY_FINFET;
							heightDrainP = maxNumPSheet * tech.widthFin; // modified from heightDrainP = (maxNumPFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
						}
					} else {
						unitWidthDrainP = 0;
						unitWidthSourceP = 0;
						heightDrainP = 0;
					}

					int NumNSheet  = (int)(ceil(widthNMOS));

					if (NumNSheet > 0) {
						if (NumNSheet <= maxNumNSheet) { /* No folding */
							unitWidthDrainN = tech.featureSize * MIN_GAP_BET_GATE_POLY_FINFET;
							unitWidthSourceN = unitWidthDrainN;
							heightDrainN =  NumNSheet * tech.widthFin;
						} else {    /* Folding */
							numFoldedNMOS = (int)(ceil(NumNSheet / maxNumNSheet));
							unitWidthDrainN = (int)ceil((double)(numFoldedNMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY_FINFET;
							unitWidthSourceN = (int)floor((double)(numFoldedNMOS+1)/2) * tech.featureSize * MIN_GAP_BET_GATE_POLY_FINFET;
							heightDrainN = maxNumNSheet * tech.widthFin; // modified from heightDrainN = (maxNumNFin-1) * tech.PitchFin + 2 * tech.widthFin/2;
						}
					} else {
						unitWidthDrainN = 0;
						unitWidthSourceN = 0;
						heightDrainN = 0;
					}

				}
			}

			switch (gateType) {
			case INV:
				if (widthPMOS > 0) {
					widthDrainP = unitWidthDrainP;
					// Folding=1: both drain and source has 1 side; folding=2: drain has 2 sides and source has 0 side... etc
					widthDrainSidewallP = widthDrainP * 2 + heightDrainP * (1+(numFoldedPMOS+1)%2);
				}
				if (widthNMOS > 0) {
					widthDrainN = unitWidthDrainN;
					widthDrainSidewallN = widthDrainN * 2 + heightDrainN * (1+(numFoldedPMOS+1)%2);
				}
				break;
			case NOR:
				// If either PMOS or NMOS has folding, there is no source/drain sharing among different PMOS and NMOS devices
				if (numFoldedPMOS == 1 && numFoldedNMOS == 1) {
					if (widthPMOS > 0) {	// No need to consider the source capacitance in series PMOS because here source and drain shares
						widthDrainP = unitWidthDrainP * numInput;
						widthDrainSidewallP = widthDrainP * 2 + heightDrainP;
					}
					if (widthNMOS > 0) {	// The number of NMOS drains is not equal to the number of NMOS in parallel because drain can share
						widthDrainN = unitWidthDrainN * (int)floor((double)(numInput+1)/2);	// Use floor: assume num of source regions >= num of drain regions
						widthDrainSidewallN = widthDrainN * 2 + heightDrainN * (1-(numInput+1)%2);
					}
				} else {
					if (widthPMOS > 0) {	// Need to consider the source capacitance in series PMOS (excluding the top one)
						widthDrainP = unitWidthDrainP * numInput + (numInput-1) * unitWidthSourceP;
						widthDrainSidewallP = widthDrainP * 2
											+ heightDrainP * (1+(numFoldedPMOS+1)%2) * numInput			// Drain sidewalls
											+ heightDrainP * (1-(numFoldedPMOS+1)%2) * (numInput-1);	// Source sidewalls
					}
					if (widthNMOS > 0) {	// Drain cannot share between different NMOS
						widthDrainN = unitWidthDrainN * numInput;
						widthDrainSidewallN = widthDrainN * 2 + heightDrainN * (1+(numFoldedNMOS+1)%2) * numInput;
					}
				}
				break;
			case NAND:
				// If either PMOS or NMOS has folding, there is no source/drain sharing among different PMOS and NMOS devices
				if (numFoldedPMOS == 1 && numFoldedNMOS == 1) {
					if (widthPMOS > 0) {  // The number of PMOS drains is not equal to the number of PMOS in parallel because drain can share
						widthDrainP = unitWidthDrainP * (int)floor((double)(numInput+1)/2); // Use floor: assume num of source regions >= num of drain regions
						widthDrainSidewallP = widthDrainP * 2 + heightDrainP * (1-(numInput+1)%2);
					}
					if (widthNMOS > 0) {  // No need to consider the source capacitance in series NMOS because here source and drain shares
						widthDrainN = unitWidthDrainN * numInput;
						widthDrainSidewallN = widthDrainN * 2 + heightDrainN;
					}
				} else {
					if (widthPMOS > 0) {  // Drain cannot share between different PMOS
						widthDrainP = unitWidthDrainP * numInput;
						widthDrainSidewallP = widthDrainP * 2 + heightDrainP * (1+(numFoldedPMOS+1)%2) * numInput;
					}
					if (widthNMOS > 0) {  // Need to consider the source capacitance in series NMOS (excluding the bottom one)
						widthDrainN = unitWidthDrainN * numInput + (numInput-1) * unitWidthSourceN;
						widthDrainSidewallN = widthDrainN * 2
											+ heightDrainN * (1+(numFoldedNMOS+1)%2) * numInput         // Drain sidewalls
											+ heightDrainN * (1-(numFoldedNMOS+1)%2) * (numInput-1);    // Source sidewalls
					}
				}
				break;
			default:
				widthDrainN = widthDrainP = widthDrainSidewallP = widthDrainSidewallN = 0;
			}
			/* Junction capacitance */
			double capDrainBottomN = widthDrainN * heightDrainN * tech.capJunction;
			double capDrainBottomP = widthDrainP * heightDrainP * tech.capJunction;

			/* Sidewall capacitance */	// FIXME
			double capDrainSidewallN, capDrainSidewallP;
			capDrainSidewallP = widthDrainSidewallP * tech.capSidewall;
			capDrainSidewallN = widthDrainSidewallN * tech.capSidewall;

			/* Drain to channel capacitance */	// FIXME
			double capDrainToChannelN = numFoldedNMOS * heightDrainN * tech.capDrainToChannel;
			double capDrainToChannelP = numFoldedPMOS * heightDrainP * tech.capDrainToChannel;
	

			// 1.4 update: junction capactiance for 14 nm and beyond
			if (tech.featureSize <= 14) {
			*(capOutput) = capDrainBottomN + capDrainBottomP;
			}

			else {
			*(capOutput) = capDrainBottomN + capDrainBottomP + capDrainSidewallN + capDrainSidewallP + capDrainToChannelN + capDrainToChannelP; 
			}
			
		}	
	
	}
}

double CalculateDrainCap(
		double width, int type,
		double heightTransistorRegion, Technology tech) {
	double drainCap = 0;
	if (type == NMOS)
		CalculateGateCapacitance(INV, 1, width, 0, heightTransistorRegion, tech, NULL, &drainCap);
	else
		CalculateGateCapacitance(INV, 1, 0, width, heightTransistorRegion, tech, NULL, &drainCap);
	return drainCap;
}

double CalculateGateLeakage(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double temperature, Technology tech) {
	int tempIndex = (int)temperature - 300;
	if ((tempIndex > 100) || (tempIndex < 0)) {
		cout<<"Error: Temperature is out of range"<<endl;
		exit(-1);
	}
	double *leakN = tech.currentOffNmos;
	double *leakP = tech.currentOffPmos;
	double leakageN, leakageP;
	switch (gateType) {
	case INV:
		leakageN = widthNMOS * leakN[tempIndex];
		leakageP = widthPMOS * leakP[tempIndex];
		return MAX(leakageN, leakageP);
	case NOR:
		leakageN = widthNMOS * leakN[tempIndex] * numInput;
		if (numInput == 2) {
			return AVG_RATIO_LEAK_2INPUT_NOR * leakageN;
		}
		else {
			return AVG_RATIO_LEAK_3INPUT_NOR * leakageN;
		}
	case NAND:
		leakageP = widthPMOS * leakP[tempIndex] * numInput;
		if (numInput == 2) {
			return AVG_RATIO_LEAK_2INPUT_NAND * leakageP;
		}
		else {
			return AVG_RATIO_LEAK_3INPUT_NAND * leakageP;
		}
	default:
		return 0.0;
	}
}

double CalculateOnResistance(double width, int type, double temperature, Technology tech) {
	double r;
	int tempIndex = (int)temperature - 300;
	if ((tempIndex > 100) || (tempIndex < 0)) {
		cout<<"Error: Temperature is out of range"<<endl;
		exit(-1);
	}
	if (type == NMOS)
		r = tech.effectiveResistanceMultiplier * tech.vdd / (tech.currentOnNmos[tempIndex] * width);
	else
		r = tech.effectiveResistanceMultiplier * tech.vdd / (tech.currentOnPmos[tempIndex] * width);
	return r;
}

double CalculateTransconductance(double width, int type, Technology tech) {
	double gm;
	double vsat;
	double widthEff;

	if (tech.featureSize >= 22){	
		if (type == NMOS) {
			vsat = MIN(tech.vdsatNmos, tech.vdd - tech.vth);
			gm = (tech.effectiveElectronMobility * tech.capOx) / 2 * width / tech.phyGateLength * vsat;
		} else {
			vsat = MIN(tech.vdsatPmos, tech.vdd - tech.vth);
			gm = (tech.effectiveHoleMobility * tech.capOx) / 2 * width / tech.phyGateLength * vsat;
		}
	} else {
		if ((tech.featureSize <= 14) && (tech.featureSize >= 3)) { // FINFET Technology
			width = width/(2 * tech.featureSize);
			widthEff = int(ceil(width))*(2*tech.heightFin + tech.widthFin);
		} else {  // GAAFET Technology
			width = width/(2 * tech.featureSize);
			widthEff = int(ceil(width))*tech.effective_width*tech.max_sheet_num/tech.max_fin_per_GAA;       
		}

		if (type == NMOS) {
			gm = tech.gm_oncurrent*widthEff;	
		} else { //type==PMOS
			gm = tech.gm_oncurrent*widthEff;
    	}
	}
	
	return gm;
}

double horowitz(double tr, double beta, double rampInput, double *rampOutput) {
	double alpha;
	alpha = 1 / rampInput / tr;
	double vs = 0.5;	/* Normalized switching voltage */
	double result = tr * sqrt(log(vs) * log(vs) + 2 * alpha * beta * (1 - vs));
	if (rampOutput)
		*rampOutput = (1 - vs) / result;
	return result;
}

double CalculateWireResistance(
		double resistivity, double wireWidth, double wireThickness,
		double barrierThickness, double dishingThickness, double alphaScatter) {
	return(alphaScatter * resistivity / (wireThickness - barrierThickness - dishingThickness)
			/ (wireWidth - 2 * barrierThickness));
}

double CalculateWireCapacitance(
		double permittivity, double wireWidth, double wireThickness, double wireSpacing,
		double ildThickness, double millerValue, double horizontalDielectric,
		double verticalDielectric, double fringeCap) {
	double verticalCap, sidewallCap;
	verticalCap = 2 * permittivity * verticalDielectric * wireWidth / ildThickness;
	sidewallCap = 2 * permittivity * millerValue * horizontalDielectric * wireThickness / wireSpacing;
	return (verticalCap + sidewallCap + fringeCap);
}

// 1.4 update: add on resistance calculation without effective resistance multiplier 

double CalculateOnResistance_normal(double width, int type, double temperature, Technology tech) {
    double r;
    int tempIndex = (int)temperature - 300;
    if ((tempIndex > 100) || (tempIndex < 0)) {
        cout<<"Error: Temperature is out of range"<<endl;
        exit(-1);
    }
	double widthEff = 0;
	if (tech.featureSize < 22  && tech.featureSize >=  3 ) { // 1.4 update : up to FinFET 7 nm
        width = width/(2 * tech.featureSize);
        widthEff = (ceil(width))*(tech.effective_width);
    } else { // 1.4 update : GAA case
        width = width/(2 * tech.featureSize);
        widthEff = (ceil(width))*tech.effective_width*tech.max_sheet_num/tech.max_fin_per_GAA;
    }

    if (type == NMOS)
        r = tech.vdd / (tech.currentOnNmos[tempIndex] * widthEff);
    else
        r =  tech.vdd / (tech.currentOnPmos[tempIndex] * widthEff);
	
    return r;
}
