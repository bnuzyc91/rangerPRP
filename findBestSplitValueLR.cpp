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

//to do construct
//in construct sampleIDs=0 is wrong
findBestSplitValueLR::findBestSplitValueLR(){};
//findBestSplitValueLR::findBestSplitValueLR(): data(0), sampleIDs(0), nodeID(0) , varID(0),  unique_timepoints(0), status_varID(0),response_timepointIDs(0), min_node_size(0) {
//}

//init
//init
void findBestSplitValueLR::init(Data* data,  std::vector<std::vector<size_t>>& sampleIDs, size_t nodeID ,size_t varID, std::vector<double>* unique_timepoints,
                                size_t status_varID,std::vector<size_t>* response_timepointIDs,size_t min_node_size,size_t* num_deaths, size_t* num_samples_at_risk)
{
     std:: cout << "print some1 in findBestSplitValueLR "<< std::endl;
    this->data=data;
    this->sampleIDs=sampleIDs;
    this->nodeID= nodeID;
    this->varID=varID;
    this->unique_timepoints=unique_timepoints;
    this->status_varID=status_varID;
    this->response_timepointIDs=response_timepointIDs;
    this->min_node_size=min_node_size;
    //    num_timepoints = unique_timepoints->size();
    this->num_deaths =num_deaths ;
    this->num_samples_at_risk = num_samples_at_risk;
    double value =data->get(nodeID, varID);
    size_t i=0;
    std:: cout << "inside findBestSplitValueLR value is " << value << "\n"<< std::endl;
        std:: cout << "sampleIDs is " << sampleIDs[nodeID].size() << "\n"<< std::endl;
    std:: cout << "unique_timepoints is " << unique_timepoints[i].size() << "\n"<< std::endl;
        std:: cout << "response_timepointIDs is " << response_timepointIDs[i].size() << "\n"<< std::endl;
        std:: cout << "status_varID is "<< status_varID << "\n"<< std::endl;
        std:: cout << "num_deaths is " << num_deaths[i] << "\n"<< std::endl;
        std:: cout << "num_samples_at_risk is " << num_samples_at_risk[i] << "\n"<< std::endl;
        std:: cout << "min_node_size is "<< min_node_size << "\n"<< std::endl;
    
         this->varID=varID;
        std::cout << typeid(varID).name() << std::endl;
        std:: cout << "inside findBestSplitValueLR varID is " << varID << "\n"<< std::endl;
}


//deconstruct todo
findBestSplitValueLR::~findBestSplitValueLR(){
    
    // Delete sampleID vector to save memory
    //sampleIDs.clear();
    //sampleIDs.shrink_to_fit();

}

void findBestSplitValueLR::printSome()
{
     std:: cout << "print some2 in findBestSplitValueLR "<< std::endl;
}

void findBestSplitValueLR::findBestSplitValueLogRank1(size_t nodeID, size_t varID, double& best_value,double& best_logrank) {
    
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
    
    computeChildDeathCounts1(nodeID, varID, possible_split_values, num_samples_right_child,
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
        
       // if (logrank > best_logrank) {
            best_value = (possible_split_values[i] + possible_split_values[i + 1]) / 2;
            //best_varID = varID;
            best_logrank = logrank;
            
            // Use smaller value if average is numerically the same as the larger value
            if (best_value == possible_split_values[i + 1]) {
                best_value = possible_split_values[i];
            }
      //  }
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

