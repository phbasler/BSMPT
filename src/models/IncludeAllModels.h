/*
 * IncludeAllModels.h
 *
 *  Created on: Mar 10, 2017
 *      Author: basler
 */

//#ifndef INCLUDEALLMODELS_H_
//#define INCLUDEALLMODELS_H_

#pragma once


#include "ClassPotentialOrigin.h"

#include <memory>

/**
 * Here are the Models of the already implemented models
 */

const int C_ModelC2HDM = 0;
const int C_ModelR2HDM = 1;
const int C_ModelRN2HDM = 2;
const int C_ModelTemplate=5;



std::unique_ptr<Class_Potential_Origin> FChoose(int choice);


//#endif /* INCLUDEALLMODELS_H_ */
