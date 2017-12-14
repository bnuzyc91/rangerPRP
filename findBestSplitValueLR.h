//
//  findBestSplitValueLR.h
//  
//
//  Created by yichen zhou on 12/10/17.
//
//

#ifndef ____findBestSplitValueLR__
#define ____findBestSplitValueLR__

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

#include "utility.h"
#include "TreeSurvival.h"
#include "Data.h"

class findBestSplitValueLR {
public:
    // trial and error
    findBestSplitValueLR();
    virtual ~findBestSplitValueLR();
   void init(Data* data,  std::vector<std::vector<size_t>>& sampleIDs, size_t nodeID ,size_t varID, std::vector<double>* unique_timepoints,
         size_t status_varID,std::vector<size_t>* response_timepointIDs,size_t min_node_size,std::size_t* num_deaths, size_t* num_samples_at_risk);
//    void init();
    void printSome();
    void findBestSplitValueLogRank1(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,double& best_logrank);
   void computeChildDeathCounts1(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
                             size_t* num_samples_right_child, size_t* delta_samples_at_risk_right_child, size_t* num_deaths_right_child,
                                 size_t num_splits);
    size_t status_varID;
    
    // Unique time points for all individuals (not only this bootstrap), sorted
    Data* data;
    size_t nodeID;
    size_t varID;
    std::vector<size_t>* response_timepointIDs;
    std::vector<double>* unique_timepoints;
    size_t num_timepoints;
    size_t  min_node_size;
    std::vector<std::vector<size_t>> sampleIDs;
    size_t* num_deaths;
    size_t* num_samples_at_risk;
    //std::vector<size_t>* num_deaths;
    //std::vector<size_t>* num_samples_at_risk;
      DISALLOW_COPY_AND_ASSIGN(findBestSplitValueLR);
};
#endif /* defined(____findBestSplitValueLR__) */
