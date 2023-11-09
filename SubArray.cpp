//Copyright (c) 2015-2016, UT-Battelle, LLC. See LICENSE file in the top-level directory
// This file contains code from NVSim, (c) 2012-2013,  Pennsylvania State University 
//and Hewlett-Packard Company. See LICENSE_NVSim file in the top-level directory.
//No part of DESTINY Project, including this file, may be copied,
//modified, propagated, or distributed except according to the terms
//contained in the LICENSE file.


#include "SubArray.h"
#include "formula.h"
#include "global.h"
#include "constant.h"
#include <math.h>

SubArray::SubArray() {
	// TODO Auto-generated constructor stub
	initialized = false;
	invalid = false;
}

SubArray::~SubArray() {
	// TODO Auto-generated destructor stub
}

void SubArray::Initialize(long long _numRow, long long _numColumn, bool _multipleRowPerSet, bool _split,
		int _muxSenseAmp, bool _internalSenseAmp, int _muxOutputLev1, int _muxOutputLev2,
		BufferDesignTarget _areaOptimizationLevel, int _num3DLevels) {
	if (initialized)
		cout << "[Subarray] Warning: Already initialized!" << endl;

	numRow = _numRow;
	numColumn = _numColumn;
	multipleRowPerSet = _multipleRowPerSet;
	split = _split;
	muxSenseAmp = _muxSenseAmp;
	muxOutputLev1 = _muxOutputLev1;
	muxOutputLev2 = _muxOutputLev2;
	internalSenseAmp = _internalSenseAmp;
	areaOptimizationLevel = _areaOptimizationLevel;
    num3DLevels = _num3DLevels;

	double maxWordlineCurrent = 0;
	double maxBitlineCurrent = 0;
	
	activityRowRead = activityRowWrite = 1/numRow;

	/* Check if the configuration is legal */
	if (inputParameter->designTarget == cache && inputParameter->cacheAccessMode != sequential_access_mode) {
		/* In these cases, each column should hold part of data in all the ways */
		if (numColumn < inputParameter->associativity) {
			invalid = true;
			initialized = true;
			return;
		}
	}

	if (cell->memCellType == DRAM || cell->memCellType == eDRAM) {
		if (muxSenseAmp > 1) {
			/* DRAM does not allow muxed bitline because of its destructive readout */
			invalid = true;
			initialized = true;
			return;
		}
	}

	if (cell->memCellType == SLCNAND) {
		if (numRow < inputParameter->flashBlockSize / inputParameter->pageSize) {
			/* SLC NAND does not have enough rows to hold the page count */
			invalid = true;
			initialized = true;
			return;
		}
		if (internalSenseAmp && muxSenseAmp < 2) {
			/* There is no way to put the sense amp */
			invalid = true;
			initialized = true;
			return;
		}
	}

	if (cell->memCellType == memristor || cell->memCellType == FBRAM) {
		if (internalSenseAmp && muxSenseAmp < 2) {
			/* There is no way to put the sense amp */
			invalid = true;
			initialized = true;
			return;
		}
	}

	if (cell->memCellType == FBRAM) {
		if (cell->resistanceOff / cell->resistanceOn < numRow / BITLINE_LEAKAGE_TOLERANCE) {
			/* bitline too long */
			invalid = true;
			initialized = true;
			return;
		}
		maxBitlineCurrent = MAX(cell->resetCurrent, cell->setCurrent) + cell->leakageCurrentAccessDevice * (numRow - 1);
	}

	if (cell->memCellType == MRAM || cell->memCellType == PCRAM || cell->memCellType == memristor) {
		if (cell->accessType == CMOS_access){
			if (tech->currentOnNmos[inputParameter->temperature - 300]
									/ tech->currentOffNmos[inputParameter->temperature - 300] < numRow / BITLINE_LEAKAGE_TOLERANCE) {
				/* bitline too long */
				invalid = true;
				initialized = true;
				return;
			}
			maxBitlineCurrent = MAX(cell->resetCurrent, cell->setCurrent) + cell->leakageCurrentAccessDevice * (numRow - 1);
		} else { //non-CMOS access
			/* Write half select problem limit the array size */
			double resetCurrent;
			if (cell->resetCurrent == 0) {
				resetCurrent = (fabs (cell->resetVoltage) - cell->voltageDropAccessDevice) / cell->resistanceOnAtResetVoltage;
			} else
				resetCurrent = cell->resetCurrent;
			int numSelectedColumnPerRow = numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
			if (cell->accessType == none_access) {
                // Based on Equation (1) in DATE2011 "Design Implications of Memristor-Based RRAM Cross-Point Structures" Xu et. al
				maxWordlineCurrent = resetCurrent * numSelectedColumnPerRow + resetCurrent * cell->resistanceOnAtResetVoltage
						/ 2 / cell->resistanceOnAtHalfResetVoltage * (numColumn - numSelectedColumnPerRow);
                maxWordlineCurrent += resetCurrent * cell->resistanceOnAtResetVoltage / 2 / cell->resistanceOnAtHalfResetVoltage
                                    * numColumn * (num3DLevels -1);
			} else { //diode or BJT
				maxWordlineCurrent = resetCurrent * numSelectedColumnPerRow + cell->leakageCurrentAccessDevice
						* (numColumn - numSelectedColumnPerRow);
                maxWordlineCurrent += cell->leakageCurrentAccessDevice * numColumn * (num3DLevels - 1);
			}
			double minWordlineDriverWidth = maxWordlineCurrent / tech->currentOnNmos[inputParameter->temperature - 300];
			if (minWordlineDriverWidth > inputParameter->maxNmosSize * tech->featureSize) {
				invalid = true;
				return;
			}
			if (cell->accessType == none_access) {
                // Based on Table 1, Row 1 in DATE2011 "Design Implications of Memristor-Based RRAM Cross-Point Structures" Xu et. al
				maxBitlineCurrent = resetCurrent + resetCurrent * cell->resistanceOnAtResetVoltage / 2
						/ cell->resistanceOnAtHalfResetVoltage * (numRow - 1);
                maxBitlineCurrent = resetCurrent * cell->resistanceOnAtResetVoltage / 2 / cell->resistanceOnAtHalfResetVoltage
                                  * numRow * (num3DLevels - 1);
			} else { //diode or BJT
				maxBitlineCurrent = resetCurrent + cell->leakageCurrentAccessDevice * (numRow - 1);
                maxBitlineCurrent += cell->leakageCurrentAccessDevice * numRow * (num3DLevels - 1);
			}
		}
	}

	double minBitlineMuxWidth = maxBitlineCurrent / tech->currentOnNmos[inputParameter->temperature - 300];
	minBitlineMuxWidth = MAX(MIN_NMOS_SIZE * tech->featureSize, minBitlineMuxWidth);
	if (minBitlineMuxWidth > inputParameter->maxNmosSize * tech->featureSize) {
		invalid = true;
		return;
	}

	if (internalSenseAmp) {
		if (cell->memCellType == SRAM || cell->memCellType == DRAM || cell->memCellType == eDRAM) {
			/* SRAM, DRAM, and eDRAM all use voltage sensing */
			voltageSense = true;
		} else if (cell->memCellType == MRAM || cell->memCellType == PCRAM || cell->memCellType == memristor || cell->memCellType == FBRAM) {
			voltageSense = cell->readMode;
		} else {/* NAND flash */
			voltageSense = true;
		}
	} else if (cell->memCellType == DRAM || cell->memCellType == eDRAM) {
		cout << "[Subarray] Error: DRAM does not support external sense amplifiers!" << endl;
		exit(-1);
	}

	//if (cell->memCellType == DRAM || cell->memCellType == eDRAM) {
	//	senseVoltage = tech->vdd / 2 * cell->capDRAMCell / (cell->capDRAMCell + capBitline);
	//	if (senseVoltage < cell->minSenseVoltage) {		/* Bitline is too long */
	//		invalid = true;
	//		initialized = true;
	//		return;
	//	}
	//} else if (cell->memCellType == SLCNAND){
	//	/* suppose the reference voltage is 0.5Vdd, the initial bitline voltage is 0.6Vdd
	//	 * if the bitline drops to 0.4Vdd, the senseamp can tell which data is stored */
	//	senseVoltage = MAX(cell->minSenseVoltage, 0.2 * tech->vdd);
	//} else {
	//	/* TO-DO: different memory technology might have different values here */
	//	senseVoltage = cell->minSenseVoltage;
	//}
	double MIN_CELL_HEIGHT = MAX_TRANSISTOR_HEIGHT;  //set real layout cell height
	double MIN_CELL_WIDTH = (MIN_GAP_BET_GATE_POLY + POLY_WIDTH) * 2;  //set real layout cell width
	double ISOLATION_REGION = MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY; // 1.4 update : new variable

	/* Add cell relaxation parameters from NeurSim 1.4 */
	switch(tech->featureSizeInNano){
		case 14:
			MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_14nm/MAX_TRANSISTOR_HEIGHT);
			ISOLATION_REGION *= (OUTER_HEIGHT_REGION_14nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));
			MIN_CELL_WIDTH  *= ((POLY_WIDTH_FINFET + MIN_GAP_BET_GATE_POLY_FINFET )/(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			break;
		case 10:
			MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_10nm /MAX_TRANSISTOR_HEIGHT);
			ISOLATION_REGION *= (OUTER_HEIGHT_REGION_10nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));
			MIN_CELL_WIDTH  *= (CPP_10nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			break;
		case 7:
			MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_7nm /MAX_TRANSISTOR_HEIGHT);
			ISOLATION_REGION *= (OUTER_HEIGHT_REGION_7nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));
			MIN_CELL_WIDTH  *= (CPP_7nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			break;
		case 5:
			MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_5nm /MAX_TRANSISTOR_HEIGHT);
			ISOLATION_REGION *= (OUTER_HEIGHT_REGION_5nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));
			MIN_CELL_WIDTH  *= (CPP_5nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			break;
		case 3:
			MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_3nm /MAX_TRANSISTOR_HEIGHT);
			ISOLATION_REGION *= (OUTER_HEIGHT_REGION_3nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));
			MIN_CELL_WIDTH  *= (CPP_3nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			break;
		case 2:
			MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_2nm /MAX_TRANSISTOR_HEIGHT);
			ISOLATION_REGION *= (OUTER_HEIGHT_REGION_2nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));
			MIN_CELL_WIDTH  *= (CPP_2nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			break;
		case 1:
			MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_1nm /MAX_TRANSISTOR_HEIGHT);
			ISOLATION_REGION *= (OUTER_HEIGHT_REGION_1nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));
			MIN_CELL_WIDTH  *= (CPP_1nm/(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
			break;
	}

	/* Derived parameters */
	numSenseAmp = numColumn / muxSenseAmp;
	if(inputParameter->relaxSRAMCell){
		lenWordline = (double)numColumn * MAX(cell->widthInFeatureSize, MIN_CELL_HEIGHT) * tech->featureSize;
		lenBitline = (double)numRow * MAX(cell->heightInFeatureSize, MIN_CELL_WIDTH) * tech->featureSize;
	} else {
		lenWordline = (double)numColumn * cell->widthInFeatureSize * tech->featureSize;
		lenBitline = (double)numRow * cell->heightInFeatureSize * tech->featureSize;
	}
	/* Add stitching overhead if necessary */
	if (cell->stitching) {
		lenWordline += ((numColumn - 1) / cell->stitching + 1) * STITCHING_OVERHEAD * devtech->featureSize;
	}
	/* Add select transistors into the length calculation */
	if (cell->memCellType == SLCNAND) {
		int pageCount = inputParameter->flashBlockSize / inputParameter->pageSize;
		/* Two select transistor including contacts have total length of 5F */
		lenBitline += (numRow / pageCount) * 5 * tech->featureSize;
	}
	/* Calculate wire resistance/capacitance */
	capWordline = lenWordline * localWire->capWirePerUnit * num3DLevels;
	resWordline = lenWordline * localWire->resWirePerUnit_M1 * num3DLevels;
	capBitline = lenBitline * localWire->capWirePerUnit * num3DLevels;
	resBitline = lenBitline * localWire->resWirePerUnit_M0 * num3DLevels;

	//cout << "capWirePerUnit: " << localWire->capWirePerUnit * 1e12 << endl;
	//cout << "resWirePerUnit_M0: " << localWire->resWirePerUnit_M0 * 1e-6 << endl;
	//cout << "resWirePerUnit_M1: " << localWire->resWirePerUnit_M1 * 1e-6<< endl;

	/* Caclulate the load resistance and capacitance for Mux Decoders */
	double capMuxLoad, resMuxLoad;
        resMuxLoad = resWordline;
        capMuxLoad = CalculateGateCap(minBitlineMuxWidth, *tech) * numColumn;
        capMuxLoad += capWordline;

		//cout << "capMuxLoad: " << capMuxLoad * 1e9 << endl;
		//cout << "resMuxLoad: " << resMuxLoad * 1e9 << endl;

	if (cell->memCellType == DRAM || cell->memCellType == eDRAM) {
		senseVoltage = devtech->vdd / 2 * cell->capDRAMCell / (cell->capDRAMCell + capBitline);
		if (senseVoltage < cell->minSenseVoltage) {		/* Bitline is too long */
			invalid = true;
			initialized = true;
			return;
		}
	} else if (cell->memCellType == SLCNAND){
		/* suppose the reference voltage is 0.5Vdd, the initial bitline voltage is 0.6Vdd
		 * if the bitline drops to 0.4Vdd, the senseamp can tell which data is stored */
		senseVoltage = MAX(cell->minSenseVoltage, 0.2 * tech->vdd);
	} else {
		/* TO-DO: different memory technology might have different values here */
		senseVoltage = cell->minSenseVoltage;
	}

	/* Add transistor resistance/capacitance */
	if (cell->memCellType == SRAM) {
		/* SRAM has two access transistors */
		resCellAccess = CalculateOnResistance(((tech->featureSize <= 14*1e-9)? 2:1) * cell->widthAccessCMOS * tech->featureSize, NMOS, inputParameter->temperature, *tech);
		capCellAccess = CalculateDrainCap(((tech->featureSize <= 14*1e-9)? 2:1) * cell->widthAccessCMOS * tech->featureSize, NMOS, cell->widthInFeatureSize * tech->featureSize, *tech);
		capWordline += 2 * CalculateGateCap(((tech->featureSize <= 14*1e-9)? 2:1) * cell->widthAccessCMOS * tech->featureSize, *tech) * numColumn;
		if(tech->featureSize <= 14 * 1e-9){ capBitline += tech->cap_draintotal * cell->widthAccessCMOS * tech->effective_width * numRow / 2;}
		else {capBitline  += capCellAccess * numRow / 2;	/* Due to shared contact */}
		voltagePrecharge = tech->vdd / 2;	/* SRAM read voltage is always half of vdd */
	} else if (cell->memCellType == DRAM || cell->memCellType == eDRAM) {
		/* DRAM and eDRAM only has one access transistors */
		resCellAccess = CalculateOnResistance(((tech->featureSize <= 14*1e-9)? 2:1) * cell->widthAccessCMOS * devtech->featureSize, NMOS, inputParameter->temperature, *devtech);
		capCellAccess = CalculateDrainCap(((tech->featureSize <= 14*1e-9)? 2:1) * cell->widthAccessCMOS * devtech->featureSize, NMOS, cell->widthInFeatureSize * devtech->featureSize, *devtech);
		capWordline += CalculateGateCap(((tech->featureSize <= 14*1e-9)? 2:1) * cell->widthAccessCMOS * devtech->featureSize, *devtech) * numColumn;
		capBitline  += capCellAccess * numRow / 2;	/* Due to shared contact */
		voltagePrecharge = devtech->vdd / 2;	/* DRAM read voltage is always half of vdd */
	} else if (cell->memCellType == FBRAM){ /* Floating Body RAM */
		resCellAccess = 0;
		capCellAccess = CalculateFBRAMDrainCap(cell->widthSOIDevice * tech->featureSize, *tech);
		capWordline += CalculateFBRAMGateCap(cell->widthSOIDevice * tech->featureSize, cell->gateOxThicknessFactor, *tech) * numColumn;
		capBitline  += capCellAccess * numRow / 2;	/* Due to shared contact */
		resMemCellOff = cell->resistanceOff;
		resMemCellOn = cell->resistanceOn;
		if (cell->readMode) {						/* voltage-sensing */
			if (cell->readVoltage == 0) {  /* Current-in voltage sensing */
				voltageMemCellOff = cell->readCurrent * resMemCellOff;
				voltageMemCellOn = cell->readCurrent * resMemCellOn;
				voltagePrecharge = (voltageMemCellOff + voltageMemCellOn) / 2;
				voltagePrecharge = MIN(tech->vdd, voltagePrecharge);  /* TO-DO: we can have charge bump to increase SA working point */
				if ((voltagePrecharge - voltageMemCellOn) <= senseVoltage) {
					cout <<"Error[Subarray]: Read current too large or too small that no reasonable precharge voltage existing" <<endl;
					invalid = true;
					return;
				}
			} else {   /*Voltage-divider sensing */
				resInSerialForSenseAmp = sqrt(resMemCellOn * resMemCellOff);
				resEquivalentOn = resMemCellOn * resInSerialForSenseAmp / (resMemCellOn + resInSerialForSenseAmp);
				resEquivalentOff = resMemCellOff * resInSerialForSenseAmp / (resMemCellOff + resInSerialForSenseAmp);
				voltageMemCellOff = cell->readVoltage * resMemCellOff / (resMemCellOff + resInSerialForSenseAmp);
				voltageMemCellOn = cell->readVoltage * resMemCellOn / (resMemCellOn + resInSerialForSenseAmp);
				voltagePrecharge = (voltageMemCellOff + voltageMemCellOn) / 2;
				voltagePrecharge = MIN(tech->vdd, voltagePrecharge);  /* TO-DO: we can have charge bump to increase SA working point */
				if ((voltagePrecharge - voltageMemCellOn) <= senseVoltage) {
					cout <<"Error[Subarray]: Read Voltage too large or too small that no reasonable precharge voltage existing" <<endl;
					invalid = true;
					return;
				}
			}
		}
	} else if (cell->memCellType == MRAM || cell->memCellType == PCRAM || cell->memCellType == memristor) {
		/* MRAM, PCRAM, and memristor have three types of access devices: CMOS, BJT, and diode */
		if (cell->accessType == CMOS_access) {
			resCellAccess = CalculateOnResistance(((tech->featureSize <= 14*1e-9)? 2:1) * cell->widthAccessCMOS * tech->featureSize, NMOS, inputParameter->temperature, *tech);
			capCellAccess = CalculateDrainCap(((tech->featureSize <= 14*1e-9)? 2:1) * cell->widthAccessCMOS * tech->featureSize, NMOS, cell->widthInFeatureSize * tech->featureSize, *tech);
			capWordline += CalculateGateCap(((tech->featureSize <= 14*1e-9)? 2:1) * cell->widthAccessCMOS * tech->featureSize, *tech) * numColumn;
			capBitline  += capCellAccess * numRow / 2;	/* Due to shared contact */
		} else if (cell->accessType == BJT_access) {
			// TO-DO
	/*	} else if (cell->accessType == diode_access){
			if (cell->readVoltage == 0) {
				resCellAccess = cell->voltageDropAccessDevice / cell->readCurrent;
			} else {
				if (cell->readMode == false) {
					resCellAccess = cell->voltageDropAccessDevice / (cell->readVoltage
							- cell->voltageDropAccessDevice) * cell->resistanceOn;
				} else {
					cout<<"Error[Subarray]: Diode access do not support voltage-input voltage sensing" <<endl;
					exit(-1);
				}
			}
			capCellAccess = MAX(cell->capacitanceOn, cell->capacitanceOff);
			capWordline += MAX(cell->capacitanceOff, cell->capacitanceOn) * numColumn;
			capBitline += MAX(cell->capacitanceOff, cell->capacitanceOn) * numRow;      */
		} else { // none_access || diode_access
			resCellAccess = 0;
			capCellAccess = MAX(cell->capacitanceOn, cell->capacitanceOff);
			capWordline += MAX(cell->capacitanceOff, cell->capacitanceOn) * numColumn;  //TO-DO: choose the right capacitance
			capBitline += MAX(cell->capacitanceOff, cell->capacitanceOn) * numRow;      //TO-DO: choose the right capacitance

            // Add capacitance for other monolithic layers
			capWordline += MAX(cell->capacitanceOff, cell->capacitanceOn) * numColumn * (num3DLevels-1);  //TO-DO: choose the right capacitance
			capBitline += MAX(cell->capacitanceOff, cell->capacitanceOn) * numRow * (num3DLevels-1);      //TO-DO: choose the right capacitance
		}
		resMemCellOff = resCellAccess + cell->resistanceOff;
		resMemCellOn = resCellAccess + cell->resistanceOn;
		if (cell->readMode) {						/* voltage-sensing */
			if (cell->readVoltage == 0) {  /* Current-in voltage sensing */
				voltageMemCellOff = cell->readCurrent * resMemCellOff;
				voltageMemCellOn = cell->readCurrent * resMemCellOn;
				voltagePrecharge = (voltageMemCellOff + voltageMemCellOn) / 2;
				voltagePrecharge = MIN(tech->vdd, voltagePrecharge);  /* TO-DO: we can have charge bump to increase SA working point */
				if ((voltagePrecharge - voltageMemCellOn) <= senseVoltage) {
					cout <<"Error[Subarray]: Read current too large or too small that no reasonable precharge voltage existing" <<endl;
					invalid = true;
					return;
				}
			} else {   /*Voltage-in voltage sensing */
				resInSerialForSenseAmp = sqrt(resMemCellOn * resMemCellOff);
				resEquivalentOn = resMemCellOn * resInSerialForSenseAmp / (resMemCellOn + resInSerialForSenseAmp);
				resEquivalentOff = resMemCellOff * resInSerialForSenseAmp / (resMemCellOff + resInSerialForSenseAmp);
				voltageMemCellOff = cell->readVoltage * resMemCellOff / (resMemCellOff + resInSerialForSenseAmp);
				voltageMemCellOn = cell->readVoltage * resMemCellOn / (resMemCellOn + resInSerialForSenseAmp);
				voltagePrecharge = (voltageMemCellOff + voltageMemCellOn) / 2;
				voltagePrecharge = MIN(tech->vdd, voltagePrecharge);  /* TO-DO: we can have charge bump to increase SA working point */
				if ((voltagePrecharge - voltageMemCellOn) <= senseVoltage) {
					cout <<"Error[Subarray]: Read Voltage too large or too small that no reasonable precharge voltage existing" <<endl;
					invalid = true;
					return;
				}
			}
		}
	} else if (cell->memCellType == SLCNAND) {
		/* Calculate the NAND flash string length, which is the page count per block plus 2 (two select transistors) */
		int pageCount = inputParameter->flashBlockSize / inputParameter->pageSize;
		int stringLength = pageCount + 2;
		resCellAccess = CalculateOnResistance(((tech->featureSize <= 14*1e-9)? 2:1) * tech->featureSize, NMOS, inputParameter->temperature, *tech) * stringLength;
		capCellAccess = CalculateDrainCap(((tech->featureSize <= 14*1e-9)? 2:1) * tech->featureSize, NMOS, cell->widthInFeatureSize * tech->featureSize, *tech);
		/* The capacitance of each cell at the gate terminal is the series of C_control_gate | C_floating_gate */
		capWordline += CalculateGateCap(((tech->featureSize <= 14*1e-9)? 2:1) * tech->featureSize, *tech) * numColumn * cell->gateCouplingRatio / (cell->gateCouplingRatio + 1);
		capBitline  += capCellAccess * (numRow / pageCount) / 2;	/* 2 is due to shared contact and the effective row count is numRow/pageCount */
		voltagePrecharge = tech->vdd * 0.6;	/* SLC NAND flash bitline precharge voltage is assumed to 0.6Vdd */
	} else {	/* MLC NAND flash */
		// TO-DO
	}

	/* Repeater Insertion Scheme */

	if (inputParameter->numRepeaters > 0) {

			sectionres = resWordline / (inputParameter->numRepeaters + 1);
			sectioncap = capWordline / (inputParameter->numRepeaters + 1);
			targetdriveres = CalculateOnResistance(((tech->featureSize <= 14*1e-9)? 2:1) * tech->featureSize * inputParameter->bufferSizeRatio, NMOS, inputParameter->temperature, *tech) ;
			widthInvN  = CalculateOnResistance(((tech->featureSize <= 14*1e-9)? 2:1) * tech->featureSize, NMOS, inputParameter->temperature, *tech) / targetdriveres * tech->featureSize;
			widthInvP = CalculateOnResistance(((tech->featureSize <= 14*1e-9)? 2:1) * tech->featureSize, PMOS, inputParameter->temperature, *tech) / targetdriveres * tech->featureSize ;
			if (tech->featureSize <= 14*1e-9){
				widthInvN = 2* ceil(widthInvN / tech->featureSize) * tech->featureSize;
				widthInvP = 2* ceil(widthInvP / tech->featureSize) * tech->featureSize;
			}

	} else {
		
		sectionres = resWordline;
		sectioncap = capWordline;

	}

	gateCapRep = CalculateGateCap(((tech->featureSize <= 14*1e-9)? 2:1) * widthInvN * tech->featureSize, *tech) + CalculateGateCap(((tech->featureSize <= 14*1e-9)? 2:1) * widthInvP * tech->featureSize, *tech);

	/****************************/

	/* Initialize sub-component */

	precharger.Initialize(tech->vdd, numColumn, capBitline, resBitline, lenBitline);
	precharger.CalculateRC();

	rowDecoder.Initialize(numRow, (sectioncap + gateCapRep), sectionres, multipleRowPerSet, areaOptimizationLevel, maxWordlineCurrent, true, lenWordline);
	if (rowDecoder.invalid) {
		invalid = true;
		return;
	}
	rowDecoder.CalculateRC();

	if (!invalid) {
		bitlineMuxDecoder.Initialize(muxSenseAmp, capMuxLoad, resMuxLoad /* TO-DO: need to fix */, false, areaOptimizationLevel, 0, true, lenWordline);
		if (bitlineMuxDecoder.invalid)
			invalid = true;
		else
			bitlineMuxDecoder.CalculateRC();
	}

	if (!invalid) { 
		senseAmpMuxLev1Decoder.Initialize(muxOutputLev1, capMuxLoad, resMuxLoad /* TO-DO: need to fix */, false, areaOptimizationLevel, 0, true, lenWordline);
		if (senseAmpMuxLev1Decoder.invalid)
			invalid = true;
		else
			senseAmpMuxLev1Decoder.CalculateRC();
	}

	if (!invalid) {
		senseAmpMuxLev2Decoder.Initialize(muxOutputLev2, capMuxLoad, resMuxLoad /* TO-DO: need to fix */, false, areaOptimizationLevel, 0, true, lenWordline);
		if (senseAmpMuxLev2Decoder.invalid)
			invalid = true;
		else
			senseAmpMuxLev2Decoder.CalculateRC();
	}

	senseAmpMuxLev2.Initialize(muxOutputLev2, numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2, 0, 0 /* TO-DO: need to fix */, maxBitlineCurrent);
	senseAmpMuxLev2.CalculateRC();

	senseAmpMuxLev1.Initialize(muxOutputLev1, numColumn / muxSenseAmp / muxOutputLev1,
			senseAmpMuxLev2.capForPreviousDelayCalculation, senseAmpMuxLev2.capForPreviousPowerCalculation, maxBitlineCurrent);
	senseAmpMuxLev1.CalculateRC();

	if (internalSenseAmp) {
		if (!invalid) {
			senseAmp.Initialize(numSenseAmp, !voltageSense, senseVoltage, lenWordline / numColumn * muxSenseAmp);
			if (senseAmp.invalid)
				invalid = true;
			else
				senseAmp.CalculateRC();
		}
		if (!invalid) {
			bitlineMux.Initialize(muxSenseAmp, numColumn / muxSenseAmp, senseAmp.capLoad, senseAmp.capLoad, maxBitlineCurrent);
		}
	} else {
		if (!invalid) {
			bitlineMux.Initialize(muxSenseAmp, numColumn / muxSenseAmp,
					senseAmpMuxLev1.capForPreviousDelayCalculation, senseAmpMuxLev1.capForPreviousPowerCalculation, maxBitlineCurrent);
		}
	}

	if (!invalid) {
		bitlineMux.CalculateRC();
	}

	initialized = true;
}

void SubArray::CalculateArea() {
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else if (invalid) {
		height = width = area = invalid_value;
	} else {
		
		double addWidth = 0, addHeight = 0;
		double bufferarea= hInv * wInv * inputParameter->numRepeaters * 2 * numRow;

		width = lenWordline;
		height = lenBitline;

		rowDecoder.CalculateArea();
		if (rowDecoder.height > height) {
			/* assume magic folding */
			addWidth = rowDecoder.area / height;
		} else {
			/* allow white space */
			addWidth = rowDecoder.width;
		}

		precharger.CalculateArea();
		if (precharger.width > width) {
			/* assume magic folding */
			addHeight = precharger.area / precharger.width;
		} else {
			/* allow white space */
			addHeight = precharger.height;
		}

		bitlineMux.CalculateArea();
		addHeight += bitlineMux.height;

		if (internalSenseAmp) {
			senseAmp.CalculateArea();
			if (senseAmp.width > width * 1.001) {
				/* should never happen */
				cout << "[ERROR] Sense Amplifier area calculation is wrong!" << endl;
			} else {
				addHeight += senseAmp.height;
			}
		}

		senseAmpMuxLev1.CalculateArea();
		addHeight += senseAmpMuxLev1.height;

		senseAmpMuxLev2.CalculateArea();
		addHeight += senseAmpMuxLev2.height;

		/* Add Repeater Cell Area */

		//addHeight += bufferarea/lenWordline;
		addWidth += bufferarea/lenBitline;

		/**************************/

		bitlineMuxDecoder.CalculateArea();
		addWidth = MAX(addWidth, bitlineMuxDecoder.width);
		senseAmpMuxLev1Decoder.CalculateArea();
		addWidth = MAX(addWidth, senseAmpMuxLev1Decoder.width);
		senseAmpMuxLev2Decoder.CalculateArea();
		addWidth = MAX(addWidth, senseAmpMuxLev2Decoder.width);

		width += addWidth;
		height += addHeight;
		area = width * height;
	}
}

void SubArray::CalculateLatency(double _rampInput) {
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else if (invalid) {
		readLatency = writeLatency = invalid_value;
	} else {

		/* Row Decoder Repeater Calculation */
		double resPullDown;
		double capLoad;
		double tr;	/* time constant */
		double gm;	/* transconductance */
		double beta;	/* for horowitz calculation */
		double rampInput = _rampInput;
		double rampOutput = _rampInput;
		double rowDecoderRepeaterLatency = 0.0;
		
		if(inputParameter->numRepeaters){
			for(int i = 0; i < inputParameter->numRepeaters; i++){
				
				if(i > 0) capLoad = sectioncap;
				else capLoad = gateCapRep + sectioncap;
				
				resPullDown = CalculateOnResistance(((tech->featureSize <= 14*1e-9)? 2:1) * widthInvN, NMOS, inputParameter->temperature, *tech);
				tr = resPullDown * capLoad + sectioncap * sectionres / 2;
				gm = CalculateTransconductance(((tech->featureSize <= 14*1e-9)? 2:1) * widthInvN, NMOS, *tech);
				beta = 1 / (resPullDown * gm);

				rowDecoderRepeaterLatency += horowitz(tr, beta, rampInput, &rampOutput);
				rampInput = rampOutput;
			}
		}
		/************************************/

		precharger.CalculateLatency(_rampInput);
		rowDecoder.CalculateLatency(_rampInput);
		bitlineMuxDecoder.CalculateLatency(_rampInput);
		senseAmpMuxLev1Decoder.CalculateLatency(_rampInput);
		senseAmpMuxLev2Decoder.CalculateLatency(_rampInput);
		columnDecoderLatency = MAX(MAX(bitlineMuxDecoder.readLatency, senseAmpMuxLev1Decoder.readLatency), senseAmpMuxLev2Decoder.readLatency);
		double decoderLatency = MAX(rowDecoder.readLatency + rowDecoderRepeaterLatency, columnDecoderLatency);
		/*need a second thought on this equation*/
		double capPassTransistor = bitlineMux.capNMOSPassTransistor +
				senseAmpMuxLev1.capNMOSPassTransistor + senseAmpMuxLev2.capNMOSPassTransistor;
		double resPassTransistor = bitlineMux.resNMOSPassTransistor +
				senseAmpMuxLev1.resNMOSPassTransistor + senseAmpMuxLev2.resNMOSPassTransistor;
		double tauChargeLatency = resPassTransistor * (capPassTransistor + capBitline) + resBitline * capBitline / 2;
		chargeLatency = horowitz(tauChargeLatency, 0, 1e20, NULL);

		if (cell->memCellType == SRAM) {
			/* Codes below calculate the bitline latency */
			double resPullDown = CalculateOnResistance(cell->widthSRAMCellNMOS * ((tech->featureSize <= 14*1e-9)? 2:1) * tech->featureSize, NMOS,
					inputParameter->temperature, *tech);
			double tau = (resCellAccess + resPullDown) * (capCellAccess + capBitline + bitlineMux.capForPreviousDelayCalculation)
					+ resBitline * (bitlineMux.capForPreviousDelayCalculation + capBitline / 2);
			tau *= log(voltagePrecharge / (voltagePrecharge - senseVoltage / 2));	/* one signal raises and the other drops, so senseVoltage/2 is enough */
			double gm = CalculateTransconductance(((tech->featureSize <= 14*1e-9)? 2:1) * cell->widthAccessCMOS * tech->featureSize, NMOS, *tech);
			double beta = 1 / (resPullDown * gm);
			double bitlineRamp = 0;
			// From NeuroSim - Elmore BL (from IDRS 2022: More Moore)
			double BLCap_perCell = capBitline / numRow + capCellAccess; // Anni update
			double BLRes_perCell = resBitline / numRow;
			double Elmore_BL = (resCellAccess + resPullDown) * BLCap_perCell * numRow   + BLCap_perCell * BLRes_perCell * numRow  * ( numRow +1 )  /2;
			bitlineDelay = Elmore_BL * log(tech->vdd / (tech->vdd - cell->minSenseVoltage / 2));
			//bitlineDelay = horowitz(tau, beta, rowDecoder.rampOutput, &bitlineRamp);
			bitlineMux.CalculateLatency(bitlineRamp);
			if (internalSenseAmp) {
				senseAmp.CalculateLatency(bitlineMuxDecoder.rampOutput);
				senseAmpMuxLev1.CalculateLatency(1e20);
				senseAmpMuxLev2.CalculateLatency(senseAmpMuxLev1.rampOutput);
			} else {
				senseAmpMuxLev1.CalculateLatency(bitlineMux.rampOutput);
				senseAmpMuxLev2.CalculateLatency(senseAmpMuxLev1.rampOutput);
			}
			readLatency = decoderLatency + bitlineDelay + bitlineMux.readLatency + senseAmp.readLatency
					+ senseAmpMuxLev1.readLatency + senseAmpMuxLev2.readLatency;
			/* assume symmetric read/write for SRAM bitline delay */
			writeLatency = readLatency;
		} else if (cell->memCellType == DRAM || cell->memCellType == eDRAM) {
			double cap = (capCellAccess + cell->capDRAMCell) * (capBitline + bitlineMux.capForPreviousDelayCalculation)
					/ (capCellAccess + cell->capDRAMCell + capBitline + bitlineMux.capForPreviousDelayCalculation);
			double res = resBitline + resCellAccess;
			double tau = 2.3 * res * cap;
			double bitlineRamp = 0;
			bitlineDelay = horowitz(tau, 0, rowDecoder.rampOutput, &bitlineRamp);
			senseAmp.CalculateLatency(bitlineRamp);
			senseAmpMuxLev1.CalculateLatency(1e20);
			senseAmpMuxLev2.CalculateLatency(senseAmpMuxLev1.rampOutput);

            /* Refresh operation does not pass sense amplifier. */
            refreshLatency = decoderLatency + bitlineDelay + senseAmp.readLatency;
            refreshLatency *= numRow; // TOTAL refresh latency for subarray
			readLatency = decoderLatency + bitlineDelay + senseAmp.readLatency
					+ senseAmpMuxLev1.readLatency + senseAmpMuxLev2.readLatency;
			/* assume symmetric read/write for DRAM/eDRAM bitline delay */
			writeLatency = readLatency;
		} else if (cell->memCellType == MRAM || cell->memCellType == PCRAM || cell->memCellType == memristor || cell->memCellType == FBRAM) {
			double bitlineRamp = 0;
			if (cell->readMode == false) {	/* current-sensing */
				/* Use ICCAD 2009 model */
				double tau = resBitline * capBitline / 2 * (resMemCellOff + resBitline / 3) / (resMemCellOff + resBitline);
                //tau *= 500.0;
				bitlineDelay = horowitz(tau, 0, rowDecoder.rampOutput, &bitlineRamp);
			} else {						/* voltage-sensing */
				if (cell->readVoltage == 0) {  /* Current-in voltage sensing */
					double tau = resMemCellOn * (capCellAccess + capBitline + bitlineMux.capForPreviousDelayCalculation)
							+ resBitline * (bitlineMux.capForPreviousDelayCalculation + capBitline / 2); /* time constant of LRS */
					bitlineDelayOn = tau * log((voltagePrecharge - voltageMemCellOn)/(voltagePrecharge - voltageMemCellOn - senseVoltage));  /* BitlineDelay of HRS */
					tau = resMemCellOff * (capCellAccess + capBitline + bitlineMux.capForPreviousDelayCalculation)
							+ resBitline * (bitlineMux.capForPreviousDelayCalculation + capBitline / 2);  /* time constant of HRS */
					bitlineDelayOff = tau * log((voltageMemCellOff - voltagePrecharge)/(voltageMemCellOff - voltagePrecharge - senseVoltage));
					bitlineDelay = MAX(bitlineDelayOn, bitlineDelayOff);
				} else {   /*Voltage-in voltage sensing */
					double tau = resEquivalentOn * (capCellAccess + capBitline + bitlineMux.capForPreviousDelayCalculation)
							+ resBitline * (bitlineMux.capForPreviousDelayCalculation + capBitline / 2); /* time constant of LRS */
					bitlineDelayOn = tau * log((voltagePrecharge - voltageMemCellOn)/(voltagePrecharge - voltageMemCellOn - senseVoltage));  /* BitlineDelay of HRS */

					tau = resEquivalentOff * (capCellAccess + capBitline + bitlineMux.capForPreviousDelayCalculation)
							+ resBitline * (bitlineMux.capForPreviousDelayCalculation + capBitline / 2);  /* time constant of HRS */
					bitlineDelayOff = tau * log((voltageMemCellOff - voltagePrecharge)/(voltageMemCellOff - voltagePrecharge - senseVoltage));
					bitlineDelay = MAX(bitlineDelayOn, bitlineDelayOff);
				}
			}
			bitlineMux.CalculateLatency(bitlineRamp);
			if (internalSenseAmp) {
				senseAmp.CalculateLatency(bitlineMuxDecoder.rampOutput);
				senseAmpMuxLev1.CalculateLatency(1e20);
				senseAmpMuxLev2.CalculateLatency(senseAmpMuxLev1.rampOutput);
			} else {
				senseAmpMuxLev1.CalculateLatency(bitlineMux.rampOutput);
				senseAmpMuxLev2.CalculateLatency(senseAmpMuxLev1.rampOutput);
			}
			readLatency = decoderLatency + bitlineDelay + bitlineMux.readLatency + senseAmp.readLatency
					+ senseAmpMuxLev1.readLatency + senseAmpMuxLev2.readLatency;

			if (cell->memCellType == PCRAM) {
				if (inputParameter->writeScheme == write_and_verify) {
					/*TO-DO: write and verify programming */
				} else {
					writeLatency = MAX(rowDecoder.writeLatency, columnDecoderLatency + chargeLatency);	/* TO-DO: why not directly use precharger latency? */
					resetLatency = writeLatency + cell->resetPulse;
					setLatency = writeLatency + cell->setPulse;
					writeLatency += MAX(cell->resetPulse, cell->setPulse);
				}
			} else if (cell->memCellType == FBRAM) {
				writeLatency = MAX(rowDecoder.writeLatency, columnDecoderLatency + chargeLatency);
				resetLatency = writeLatency + cell->resetPulse;
				setLatency = writeLatency + cell->setPulse;
				writeLatency += MAX(cell->resetPulse, cell->setPulse);
			} else { //memristor and MRAM
				if (cell->accessType == diode_access || cell->accessType == none_access) {
					if (inputParameter->writeScheme == erase_before_reset || inputParameter->writeScheme == erase_before_set)
						writeLatency = MAX(rowDecoder.writeLatency, chargeLatency);
					else
						writeLatency = MAX(rowDecoder.writeLatency, columnDecoderLatency + chargeLatency);
					writeLatency += chargeLatency;
					writeLatency += cell->resetPulse + cell->setPulse;
				} else { // CMOS or Bipolar access
					writeLatency = MAX(rowDecoder.writeLatency, columnDecoderLatency + chargeLatency);
					resetLatency = writeLatency + cell->resetPulse;
					setLatency = writeLatency + cell->setPulse;
					writeLatency += MAX(cell->resetPulse, cell->setPulse);
				}
			}
		} else if (cell->memCellType == SLCNAND) {
			/* Calculate the NAND flash string length, which is the page count per block plus 2 (two select transistors) */
			int pageCount = inputParameter->flashBlockSize / inputParameter->pageSize;
			int stringLength = pageCount + 2;
			/* Codes below calculate the bitline latency */
			double resPullDown = CalculateOnResistance(((tech->featureSize <= 14*1e-9)? 2:1) * tech->featureSize, NMOS, inputParameter->temperature, *tech)
					* stringLength;
			double tau = resPullDown * (capCellAccess + capBitline + bitlineMux.capForPreviousDelayCalculation)
					+ resBitline * (bitlineMux.capForPreviousDelayCalculation + capBitline / 2);
			/* in one case the bitline is unchanged, and in the other case the bitline drops from 0.6V to 0.4V */
			tau *= log((voltagePrecharge)/ (voltagePrecharge - senseVoltage));
			double gm = CalculateTransconductance(((tech->featureSize <= 14*1e-9)? 2:1) * tech->featureSize, NMOS, *tech);	/* minimum size transistor */
			double beta = 1 / (resPullDown * gm);
			double bitlineRamp = 0;
			bitlineDelay = horowitz(tau, beta, rowDecoder.rampOutput, &bitlineRamp);
			/* to correct unnecessary horowitz calculation, TO-DO: need to revisit */
			bitlineDelay = MAX(bitlineDelay, tau * 20);
			bitlineMux.CalculateLatency(bitlineRamp);
			if (internalSenseAmp) {
				senseAmp.CalculateLatency(bitlineMuxDecoder.rampOutput);
				senseAmpMuxLev1.CalculateLatency(1e20);
				senseAmpMuxLev2.CalculateLatency(senseAmpMuxLev1.rampOutput);
			} else {
				senseAmpMuxLev1.CalculateLatency(bitlineMux.rampOutput);
				senseAmpMuxLev2.CalculateLatency(senseAmpMuxLev1.rampOutput);
			}
			readLatency = decoderLatency + bitlineDelay + bitlineMux.readLatency + senseAmp.readLatency
					+ senseAmpMuxLev1.readLatency + senseAmpMuxLev2.readLatency;
			/* calculate the erase time, a.k.a. reset here */
			resetLatency = MAX(rowDecoder.readLatency, columnDecoderLatency + chargeLatency) + cell->flashEraseTime;
			/* calculate the programming time, a.k.a. set here */
			setLatency = MAX(rowDecoder.readLatency, columnDecoderLatency + chargeLatency) + cell->flashProgramTime;
			/* use the programming latency as the write latency for SLC NAND*/
			writeLatency = setLatency;
		} else {	/* MLC NAND */
			/* TO-DO */
		}
	}
}

void SubArray::CalculatePower() {
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else if (invalid) {
		readDynamicEnergy = writeDynamicEnergy = leakage = invalid_value;
	} else {
		precharger.CalculatePower();
		rowDecoder.CalculatePower();
		bitlineMuxDecoder.CalculatePower();
		senseAmpMuxLev1Decoder.CalculatePower();
		senseAmpMuxLev2Decoder.CalculatePower();
		bitlineMux.CalculatePower();
		if (internalSenseAmp) {
			senseAmp.CalculatePower();
		}
		senseAmpMuxLev1.CalculatePower();
		senseAmpMuxLev2.CalculatePower();

		double repeater_leakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter->temperature, *tech) * inputParameter->numRepeaters *2 * numRow;
		double repeater_readDynamicEnergy = (drivecapin + drivecapout) * tech->vdd * tech->vdd * inputParameter->numRepeaters * 2 * numRow * activityRowRead;
		double repeater_writeDynamicEnergy = (drivecapin + drivecapout) * tech->vdd * tech->vdd * inputParameter->numRepeaters * 2 * numRow * activityRowWrite;

		if (cell->memCellType == SRAM) {
			/* Codes below calculate the SRAM bitline power */
			readDynamicEnergy = (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation)
					* voltagePrecharge * voltagePrecharge * numColumn;
			writeDynamicEnergy = (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation)
					* voltagePrecharge * voltagePrecharge * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
			leakage = CalculateGateLeakage(INV, 1, cell->widthSRAMCellNMOS * tech->featureSize,
					cell->widthSRAMCellPMOS * tech->featureSize, inputParameter->temperature, *tech)
					* tech->vdd * 2;	/* two inverters per SRAM cell */
			leakage += CalculateGateLeakage(INV, 1, cell->widthAccessCMOS * tech->featureSize, 0,
					inputParameter->temperature, *tech) * tech->vdd;	/* two accesses NMOS, but combined as one with vdd crossed */
			leakage *= numRow * numColumn;

			/* Add Energy Comsumption from Repeaters */
			leakage += repeater_leakage;
			readDynamicEnergy += repeater_readDynamicEnergy;
			writeDynamicEnergy += repeater_writeDynamicEnergy;
			/*****************************************/

		} else if (cell->memCellType == DRAM || cell->memCellType == eDRAM) {
			/* Codes below calculate the DRAM bitline power */
			readDynamicEnergy = (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation) * senseVoltage * devtech->vdd * numColumn;
            refreshDynamicEnergy = readDynamicEnergy;
			double writeVoltage = cell->resetVoltage;	/* should also equal to setVoltage, for DRAM, it is Vdd */
			writeDynamicEnergy = (capBitline + bitlineMux.capForPreviousPowerCalculation) * writeVoltage * writeVoltage * numColumn;
			leakage = readDynamicEnergy / DRAM_REFRESH_PERIOD * numRow;
		} else if (cell->memCellType == MRAM || cell->memCellType == PCRAM || cell->memCellType == memristor || cell->memCellType == FBRAM) {
			if (cell->readMode == false) {	/* current-sensing */
				/* Use ICCAD 2009 model */
				double resBitlineMux = bitlineMux.resNMOSPassTransistor;
				double vpreMin = cell->readVoltage * resBitlineMux / (resBitlineMux + resBitline +resMemCellOn);
				double vpreMax = cell->readVoltage * (resBitlineMux + resBitline) / (resBitlineMux + resBitline + resMemCellOn);
				readDynamicEnergy = capCellAccess * vpreMax * vpreMax + bitlineMux.capForPreviousPowerCalculation
						* vpreMin * vpreMin + capBitline * (vpreMax * vpreMax + vpreMin * vpreMin + vpreMax * vpreMin) / 3;
				readDynamicEnergy *= numColumn;
			} else {						/* voltage-sensing */
				readDynamicEnergy = (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation) *
						(voltagePrecharge * voltagePrecharge - voltageMemCellOn * voltageMemCellOn ) * numColumn;
			}

			if (cell->readPower == 0) 
				cellReadEnergy = 2 * cell->CalculateReadPower() * senseAmp.readLatency; /* x2 is because of the reference cell */
			else
				cellReadEnergy = 2 * cell->readPower * senseAmp.readLatency;
			cellReadEnergy *= numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;

			/* Ignore the dynamic transition during the SET/RESET operation */
			/* Assume that the cell resistance keeps high for worst-case power estimation */
			cell->CalculateWriteEnergy();

			double resetEnergyPerBit = cell->resetEnergy;
			double setEnergyPerBit = cell->setEnergy;
			if (cell->setMode)
				setEnergyPerBit += (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation) * cell->setVoltage * cell->setVoltage;
			else
				setEnergyPerBit += (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation) * tech->vdd * tech->vdd;
			if (cell->resetMode)
				resetEnergyPerBit += (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation) * cell->resetVoltage * cell->resetVoltage;
			else
				resetEnergyPerBit += (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation) * tech->vdd * tech->vdd;

			if (cell->memCellType == PCRAM) { //PCRAM write energy
				if (inputParameter->writeScheme == write_and_verify) {
					/*TO-DO: write and verify programming */
				} else {
					cellResetEnergy = resetEnergyPerBit * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
					cellSetEnergy = setEnergyPerBit * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
					cellResetEnergy /= SHAPER_EFFICIENCY_CONSERVATIVE;
					cellSetEnergy /= SHAPER_EFFICIENCY_CONSERVATIVE;  /* Due to the shaper inefficiency */
					writeDynamicEnergy = MAX(cellResetEnergy, cellSetEnergy);
				}
			} else if (cell->memCellType == FBRAM){ //FBRAM write energy
				cellResetEnergy = resetEnergyPerBit * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
				cellSetEnergy = setEnergyPerBit * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
				cellResetEnergy /= SHAPER_EFFICIENCY_AGGRESSIVE;
				cellSetEnergy /= SHAPER_EFFICIENCY_AGGRESSIVE;  /* Due to the shaper inefficiency */
				writeDynamicEnergy = MAX(cellResetEnergy, cellSetEnergy);
			} else { //MRAM and memristor write energy
				if (cell->accessType == diode_access || cell->accessType == none_access) {
					if (inputParameter->writeScheme == erase_before_reset || inputParameter->writeScheme == erase_before_set) {
						cellResetEnergy = resetEnergyPerBit * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
						cellSetEnergy = setEnergyPerBit * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
						writeDynamicEnergy = cellResetEnergy + cellSetEnergy;	/* TO-DO: bug here, did you consider the write pattern? */
					} else { /* write scheme = set_before_reset or reset_before_set */
						cellResetEnergy = resetEnergyPerBit * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
						cellSetEnergy = setEnergyPerBit * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
						writeDynamicEnergy = MAX(cellResetEnergy, cellSetEnergy);
					}
				} else {
					cellResetEnergy = resetEnergyPerBit * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
					cellSetEnergy = setEnergyPerBit * numColumn / muxSenseAmp / muxOutputLev1 / muxOutputLev2;
					writeDynamicEnergy = MAX(cellResetEnergy, cellSetEnergy);
				}
				cellResetEnergy /= SHAPER_EFFICIENCY_AGGRESSIVE;
				cellSetEnergy /= SHAPER_EFFICIENCY_AGGRESSIVE;  /* Due to the shaper inefficiency */
				writeDynamicEnergy /= SHAPER_EFFICIENCY_AGGRESSIVE;
			}
			leakage = 0;                       //TO-DO: cell leaks during read/write operation
		} else if (cell->memCellType == SLCNAND) {
			/* Calculate the NAND flash string length, which is the page count per block plus 2 (two select transistors) */
			int pageCount = inputParameter->flashBlockSize / inputParameter->pageSize;
			int stringLength = pageCount + 2;

			/* === READ energy === */
			/* only the selected bitline is charged during the read operation, bitline is charged to Vpre */
			readDynamicEnergy = (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation)
					* voltagePrecharge * voltagePrecharge * numColumn;
			/* tricky thing here!
			 * In SLC NAND operation, SSL, GSL, and unselected wordlines in a block are charged to Vpass,
			 * but the selected wordline is not charged, which is totally different from the other cases.
			 */
			rowDecoder.resetDynamicEnergy = rowDecoder.readDynamicEnergy;
			rowDecoder.setDynamicEnergy = rowDecoder.readDynamicEnergy;
			double actualWordlineReadEnergy = rowDecoder.readDynamicEnergy / tech->vdd / tech->vdd
					* cell->flashPassVoltage * cell->flashPassVoltage;	/* approximate calculate, the wordline is charged to Vpass instead of Vdd */
			actualWordlineReadEnergy = actualWordlineReadEnergy * (numRow / pageCount * stringLength - 1);	/* except the selected wordline itself */
			rowDecoder.readDynamicEnergy = actualWordlineReadEnergy;	/* update the correct value */

			/* === Programming (SET) energy === */
			/* first calculate the source line energy (charged to Vdd), which is a part of "bitline" in this scenario */
			setDynamicEnergy = (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation)
					* cell->flashProgramVoltage * cell->flashProgramVoltage * numColumn;
			/* add tunneling current */
			/* originally it should be multiplied by numColumn/muxSenseAmp/muxOutputLev1/muxOutputLev2,
			 * but it is multiplied by numColumn here because all the unselected bitlines also need to precharge to Vdd
			 */
			setDynamicEnergy += DELTA_V_TH * TUNNEL_CURRENT_FLOW * cell->area
					* tech->featureSize * tech->featureSize * cell->flashProgramTime * numColumn;
			/* in programming, the SSL is precharged to Vdd, which is equal to the original value calculated
			 * from row decoder
			 */
			double actualWordlineSetEnergy = rowDecoder.setDynamicEnergy;
			/* however, the unselected wordlines in the same block have to precharge to Vpass */
			actualWordlineSetEnergy += rowDecoder.setDynamicEnergy / tech->vdd / tech->vdd
					* cell->flashPassVoltage * cell->flashPassVoltage * (numRow / pageCount * stringLength - 1);
			/* And the selected wordline is precharged to Vpgm */
			actualWordlineSetEnergy += rowDecoder.setDynamicEnergy / tech->vdd / tech->vdd
					* cell->flashProgramVoltage * cell->flashProgramVoltage;
			rowDecoder.setDynamicEnergy = actualWordlineSetEnergy;	/* update the correct value */

			/* === Erase (RESET) energy === */
			/* in erase, all the bitlines (selected or unselected) and the sourceline are precharged to Vera-Vbi */

			resetDynamicEnergy = (capCellAccess + capBitline + bitlineMux.capForPreviousPowerCalculation)
					* (cell->flashEraseVoltage - tech->buildInPotential) * (cell->flashEraseVoltage - tech->buildInPotential);
			resetDynamicEnergy *= (numColumn + 1);	/* plus 1 is due to the source line */
			/* the p-well shared by the selected block is precharged to Vera */
			double wellJunctionCap = tech->capJunction * cell->area * tech->featureSize * tech->featureSize;
			wellJunctionCap *= inputParameter->flashBlockSize;	/* one block shares the same well */
			resetDynamicEnergy += wellJunctionCap * cell->flashEraseVoltage * cell->flashEraseVoltage;
			/* in erase, all the wordlines, SSL, and GSL in unselected block are precharged to Vera * beta
			 * in selected block, SSL and GSL are precharged to Vera * beta
			 * here beta is fixed at 0.8
			 */
			double beta = 0.8;
			double actualWordlineResetEnergy = rowDecoder.resetDynamicEnergy / tech->vdd / tech->vdd
					* (cell->flashEraseVoltage * beta) * (cell->flashEraseVoltage * beta);
			actualWordlineResetEnergy *= (numRow / pageCount * stringLength - pageCount);
			rowDecoder.resetDynamicEnergy = actualWordlineResetEnergy;

			/* let write energy to be the average energy per page*/
			rowDecoder.writeDynamicEnergy = (rowDecoder.setDynamicEnergy + rowDecoder.resetDynamicEnergy / pageCount) / 2;
			writeDynamicEnergy = (setDynamicEnergy + resetDynamicEnergy / pageCount) / 2;

			/* Assume NAND flash cell does not consume any leakage */
			leakage = 0;
		} else {	/* MLC NAND */
			/* TO-DO */
		}

		if (inputParameter->designTarget == cache && inputParameter->cacheAccessMode != sequential_access_mode) {
			cellResetEnergy /= inputParameter->associativity;
			cellSetEnergy /= inputParameter->associativity;
			writeDynamicEnergy /= inputParameter->associativity;
			resetDynamicEnergy /= inputParameter->associativity;
			setDynamicEnergy /= inputParameter->associativity;
		}

		readDynamicEnergy += cellReadEnergy + rowDecoder.readDynamicEnergy + bitlineMuxDecoder.readDynamicEnergy + senseAmpMuxLev1Decoder.readDynamicEnergy
				+ senseAmpMuxLev2Decoder.readDynamicEnergy + precharger.readDynamicEnergy + bitlineMux.readDynamicEnergy
				+ senseAmp.readDynamicEnergy + senseAmpMuxLev1.readDynamicEnergy + senseAmpMuxLev2.readDynamicEnergy;
		writeDynamicEnergy += rowDecoder.writeDynamicEnergy + bitlineMuxDecoder.writeDynamicEnergy + senseAmpMuxLev1Decoder.writeDynamicEnergy
				+ senseAmpMuxLev2Decoder.writeDynamicEnergy + bitlineMux.writeDynamicEnergy
				+ senseAmp.writeDynamicEnergy + senseAmpMuxLev1.writeDynamicEnergy + senseAmpMuxLev2.writeDynamicEnergy;

        /* Read all column energy + row decoder + sense amp + precharger is enough for one subarray row REF. */
        refreshDynamicEnergy += rowDecoder.readDynamicEnergy + precharger.readDynamicEnergy
                             + senseAmp.readDynamicEnergy;
        refreshDynamicEnergy *= numRow; // Energy for this entire subarray 

		/* for assymetric RESET and SET latency calculation only */
		setDynamicEnergy += cellSetEnergy + rowDecoder.setDynamicEnergy + bitlineMuxDecoder.writeDynamicEnergy + senseAmpMuxLev1Decoder.writeDynamicEnergy
				+ senseAmpMuxLev2Decoder.writeDynamicEnergy + bitlineMux.writeDynamicEnergy
				+ senseAmp.writeDynamicEnergy + senseAmpMuxLev1.writeDynamicEnergy + senseAmpMuxLev2.writeDynamicEnergy;
		resetDynamicEnergy += setDynamicEnergy + rowDecoder.resetDynamicEnergy + bitlineMuxDecoder.writeDynamicEnergy + senseAmpMuxLev1Decoder.writeDynamicEnergy
				+ senseAmpMuxLev2Decoder.writeDynamicEnergy + bitlineMux.writeDynamicEnergy
				+ senseAmp.writeDynamicEnergy + senseAmpMuxLev1.writeDynamicEnergy + senseAmpMuxLev2.writeDynamicEnergy;

		if (cell->accessType == diode_access || cell->accessType == none_access) {
			writeDynamicEnergy += bitlineMux.writeDynamicEnergy + senseAmp.writeDynamicEnergy
					+ senseAmpMuxLev1.writeDynamicEnergy + senseAmpMuxLev2.writeDynamicEnergy;
		}
		leakage += rowDecoder.leakage + bitlineMuxDecoder.leakage + senseAmpMuxLev1Decoder.leakage
				+ senseAmpMuxLev2Decoder.leakage + precharger.leakage + bitlineMux.leakage
				+ senseAmp.leakage + senseAmpMuxLev1.leakage + senseAmpMuxLev2.leakage;
	}
}

void SubArray::PrintProperty() {
	cout << "Subarray Properties:" << endl;
	FunctionUnit::PrintProperty();
	cout << "numRow:" << numRow << " numColumn:" << numColumn << endl;
	cout << "lenWordline * lenBitline = " << lenWordline*1e6 << "um * " << lenBitline*1e6 << "um = " << lenWordline * lenBitline * 1e6 << "mm^2" << endl;
	cout << "Row Decoder Area:" << rowDecoder.height*1e6 << "um x " << rowDecoder.width*1e6 << "um = " << rowDecoder.area*1e6 << "mm^2" << endl;
	cout << "Sense Amplifier Area:" << senseAmp.height*1e6 << "um x " << senseAmp.width*1e6 << "um = " << senseAmp.area*1e6 << "mm^2" << endl;
	cout << "Subarray Area Efficiency = " << lenWordline * lenBitline / area * 100 <<"%" << endl;
	cout << "bitlineDelay: " << bitlineDelay*1e12 << "ps" << endl;
	cout << "chargeLatency: " << chargeLatency*1e12 << "ps" << endl;
	cout << "columnDecoderLatency: " << columnDecoderLatency*1e12 << "ps" << endl;
}

SubArray & SubArray::operator=(const SubArray &rhs) {
	height = rhs.height;
	width = rhs.width;
	area = rhs.area;
	readLatency = rhs.readLatency;
	writeLatency = rhs.writeLatency;
	readDynamicEnergy = rhs.readDynamicEnergy;
	writeDynamicEnergy = rhs.writeDynamicEnergy;
	resetLatency = rhs.resetLatency;
	setLatency = rhs.setLatency;
    refreshLatency = rhs.refreshLatency;
	resetDynamicEnergy = rhs.resetDynamicEnergy;
	setDynamicEnergy = rhs.setDynamicEnergy;
    refreshDynamicEnergy = rhs.refreshDynamicEnergy;
	cellReadEnergy = rhs.cellReadEnergy;
	cellResetEnergy = rhs.cellResetEnergy;
	cellSetEnergy = rhs.cellSetEnergy;
	leakage = rhs.leakage;
	initialized = rhs.initialized;
	numRow = rhs.numRow;
	numColumn = rhs.numColumn;
	multipleRowPerSet = rhs.multipleRowPerSet;
	split = rhs.split;
	muxSenseAmp = rhs.muxSenseAmp;
	internalSenseAmp = rhs.internalSenseAmp;
	muxOutputLev1 = rhs.muxOutputLev1;
	muxOutputLev2 = rhs.muxOutputLev2;
	areaOptimizationLevel = rhs.areaOptimizationLevel;
    num3DLevels = rhs.num3DLevels;

	voltageSense = rhs.voltageSense;
	senseVoltage = rhs.senseVoltage;
	numSenseAmp = rhs.numSenseAmp;
	lenWordline = rhs.lenWordline;
	lenBitline = rhs.lenBitline;
	capWordline = rhs.capWordline;
	capBitline = rhs.capBitline;
	resWordline = rhs.resWordline;
	resBitline = rhs.resBitline;
	resCellAccess = rhs.resCellAccess;
	capCellAccess = rhs.capCellAccess;
	bitlineDelay = rhs.bitlineDelay;
	chargeLatency = rhs.chargeLatency;
	columnDecoderLatency = rhs.columnDecoderLatency;
	bitlineDelayOn = rhs.bitlineDelayOn;
	bitlineDelayOff = rhs.bitlineDelayOff;
	resInSerialForSenseAmp = rhs.resInSerialForSenseAmp;
	resEquivalentOn = rhs.resEquivalentOn;
	resEquivalentOff = rhs.resEquivalentOff;
	resMemCellOff = rhs.resMemCellOff;
	resMemCellOn = rhs.resMemCellOn;

	rowDecoder = rhs.rowDecoder;
	bitlineMuxDecoder = rhs.bitlineMuxDecoder;
	bitlineMux = rhs.bitlineMux;
	senseAmpMuxLev1Decoder = rhs.senseAmpMuxLev1Decoder;
	senseAmpMuxLev1 = rhs.senseAmpMuxLev1;
	senseAmpMuxLev2Decoder = rhs.senseAmpMuxLev2Decoder;
	senseAmpMuxLev2 = rhs.senseAmpMuxLev2;
	precharger = rhs.precharger;
	senseAmp = rhs.senseAmp;

	return *this;
}
