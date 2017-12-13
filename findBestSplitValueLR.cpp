//
//  findBestSplitValueLR.cpp
//  
//
//  Created by yichen zhou on 12/10/17.
//
//
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

#include "utility.h"
#include "TreeSurvival.h"
#include "Data.h"
#include "findBestSplitValueLR.h"

//construct
findBestSplitValueLR::findBestSplitValueLR( Data* data,  std::vector<std::vector<size_t>> sampleIDs, size_t nodeID ,size_t varID, std::vector<double>* unique_timepoints,
                                           size_t status_varID,std::vector<size_t>* response_timepointIDs){
}

//init
findBestSplitValueLR::init(Data* data,  std::vector<std::vector<size_t>> sampleIDs, size_t nodeID ,size_t varID, std::vector<double>* unique_timepoints,
                           size_t status_varID,std::vector<size_t>* response_timepointIDs)
{
    this->data=data;
    this->sampleIDS=sampleIDS;
    this->nodeID= nodeID;
    this->varID=varID;
    this->unique_timepoints=unique_timepoints;
    this->status_varID=status_varID;
    this->response_timepointIDs=response_timepointIDs;
    
    
}

//deconstruct
findBestSplitValueLR::~findBestSplitValueLR(){
}

void findBestSplitValueLR::printSome()
{
     std:: cout << "print some in findBestSplitValueLR "<< std::endl;
}

void findBestSplitValueLR::findBestSplitValueLogRank1(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,double& best_logrank) {
    
    // Create possible split values
    std::vector<double> possible_split_values;
    data->getAllValues(possible_split_values, sampleIDs[nodeID], varID);
    
    // Try next variable if all equal for this
    if (possible_split_values.size() < 2) {
        return;
    }
    
    // -1 because no split possible at largest value
    size_t num_splits = possible_split_values.size() - 1;
    
    // Initialize
    size_t* num_deaths_right_child = new size_t[num_splits * num_timepoints]();
    size_t* delta_samples_at_risk_right_child = new size_t[num_splits * num_timepoints]();
    size_t* num_samples_right_child = new size_t[num_splits]();
    
    computeChildDeathCounts(nodeID, varID, possible_split_values, num_samples_right_child,
                            delta_samples_at_risk_right_child, num_deaths_right_child, num_splits);
    
    // Compute logrank test for all splits and use best
    for (size_t i = 0; i < num_splits; ++i) {
        double numerator = 0;
        double denominator_squared = 0;
        
        // Stop if minimal node size reached
        size_t num_samples_left_child = sampleIDs[nodeID].size() - num_samples_right_child[i];
        if (num_samples_right_child[i] < min_node_size || num_samples_left_child < min_node_size) {
            continue;
        }
        
        // Compute logrank test statistic for this split
        size_t num_samples_at_risk_right_child = num_samples_right_child[i];
        for (size_t t = 0; t < num_timepoints; ++t) {
            if (num_samples_at_risk[t] < 2 || num_samples_at_risk_right_child < 1) {
                break;
            }
            
            if (num_deaths[t] > 0) {
                // Numerator and demoninator for log-rank test, notation from Ishwaran et al.
                double di = (double) num_deaths[t];
                double di1 = (double) num_deaths_right_child[i * num_timepoints + t];
                double Yi = (double) num_samples_at_risk[t];
                double Yi1 = (double) num_samples_at_risk_right_child;
                numerator += di1 - Yi1 * (di / Yi);
                denominator_squared += (Yi1 / Yi) * (1.0 - Yi1 / Yi) * ((Yi - di) / (Yi - 1)) * di;
            }
            
            // Reduce number of samples at risk for next timepoint
            num_samples_at_risk_right_child -= delta_samples_at_risk_right_child[i * num_timepoints + t];
            
        }
        double logrank = -1;
        if (denominator_squared != 0) {
            logrank = fabs(numerator / sqrt(denominator_squared));
        }
        
        if (logrank > best_logrank) {
            best_value = (possible_split_values[i] + possible_split_values[i + 1]) / 2;
            best_varID = varID;
            best_logrank = logrank;
            
            // Use smaller value if average is numerically the same as the larger value
            if (best_value == possible_split_values[i + 1]) {
                best_value = possible_split_values[i];
            }
        }
    }
    
    delete[] num_deaths_right_child;
    delete[] delta_samples_at_risk_right_child;
    delete[] num_samples_right_child;
}

void findBestSplitValueLR::computeChildDeathCounts1(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
                                           size_t* num_samples_right_child, size_t* delta_samples_at_risk_right_child, size_t* num_deaths_right_child,
                                           size_t num_splits) {
    
    // Count deaths in right child per timepoint and possbile split
    for (auto& sampleID : sampleIDs[nodeID]) {
        double value = data->get(sampleID, varID);
        size_t survival_timeID = (*response_timepointIDs)[sampleID];
        
        // Count deaths until split_value reached
        for (size_t i = 0; i < num_splits; ++i) {
            
            if (value > possible_split_values[i]) {
                ++num_samples_right_child[i];
                ++delta_samples_at_risk_right_child[i * num_timepoints + survival_timeID];
                if (data->get(sampleID, status_varID) == 1) {
                    ++num_deaths_right_child[i * num_timepoints + survival_timeID];
                }
            } else {
                break;
            }
        }
    }
}

